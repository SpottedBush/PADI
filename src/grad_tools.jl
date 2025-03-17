#
# tgrad_tools.jl
#
# Provide tools for the calculus of the PADI data fidelity term. 
#
# ------------------------------------------------
#
# This file is part of PADI
#
#
# Copyright (c) 2024-2025 Vincent TARDIEUX and Laurence DENNEULIN (see LICENCE.md)
#

#------------------------------------------------

const ImageInterpolator{T<:AbstractFloat, K<:Kernel{T}} = TwoDimensionalTransformInterpolator{T,K,K}
const MyKer = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)
#const MyKer = LinearInterpolators.RectangularSpline(Float64, LinearInterpolators.Flat)
#const MyKer = LinearInterpolators.LinearSpline(Float64, LinearInterpolators.Flat)

# TODO : User should be able to select which kernel

struct parameters_table{T<: AbstractFloat}
    cols::NTuple{3, Int64} #Input size (reconstruction)
    rows::NTuple{2, Int64} #Output size (data)
    dataset_length::Int64
    Nframe::Int64
    Nrot::Int64
    Nangle::Int64
    v::Vector{NTuple{2, NTuple{3, T}}}
    indices::Array{Int64, 2}
    center::Array{T, 1}
    psf_center::NTuple{2,Array{T, 1}}
    epsilon::Vector{NTuple{2,Array{T, 1}}} 
    derotang::Vector{T}
end

struct FieldTransformOperator{T<:AbstractFloat, L<:Mapping, R<:Mapping} <: LinearMapping
    cols::NTuple{3, Int}
    rows::NTuple{2, Int}
    v_l::NTuple{3, T}
    v_r::NTuple{3, T}
    H_l_star::L
    H_l_disk::L
    H_r_star::R
    H_r_disk::R
end

struct data_table{T<:AbstractFloat, A<:FieldTransformOperator{T}} # Other types
    data::Array{T, 2};
    weights::Array{T, 2};
    H::A
end

function Set_Vi(Indices::AbstractArray{Int64,2}; alpha=[0, pi/4, pi/8, 3*pi/8], psi=[0,pi/2])
	J(a)=[cos(2*a) sin(2*a); sin(2*a) -cos(2*a)]
	P(a)=[cos(a), -sin(a)]
    v=Vector{NTuple{2, NTuple{3, Float64}}}(undef, length(Indices));
	for i=1:4
		v1=J(alpha[i])*P(psi[1])
		vn1=norm(v1)^2;
		v2=J(alpha[i])*P(psi[2])
		vn2=norm(v2)^2;
		for k=1:length(Indices[i,:])		
		v[Indices[i,k]] =((round(vn1/2, digits=4), round((abs(v1[1])^2-abs(v1[2])^2)/2, digits=4), round(real(v1[1]*v1[2]), digits=4)),
		                  (round(vn2/2, digits=4), round((abs(v2[1])^2-abs(v2[2])^2)/2, digits=4), round(real(v2[1]*v2[2]), digits=4)));
        end
	end
	return v
end

function Set_Vi(Nframe::Int, dataset_length::Int, mueller_instru::AbstractArray{Float64,2})
    @assert dataset_length == Nframe * size(mueller_instru)[1]
    
    v=Vector{NTuple{2, NTuple{3, Float64}}}(undef, dataset_length);
    for k=1:size(mueller_instru)[1]
        for l=1:Nframe
            v[Nframe*(k-1) + l] =((mueller_instru[k,1], mueller_instru[k,2], mueller_instru[k,3]),
                                  (mueller_instru[k,4], mueller_instru[k,5], mueller_instru[k,6]));
        end
    end
	return v
end


function reset_instrument(V::Vector{NTuple{2, NTuple{3, T}}})  where {T <: AbstractFloat}
    push!(Parameters, parameters_table(get_par().cols, 
                                       get_par().rows, 
                                       get_par().dataset_length, 
                                       get_par().Nframe, 
                                       get_par().Nrot, 
                                       get_par().Nangle, 
                                       V, 
                                       get_par().indices, 
                                       get_par().center, 
                                       get_par().psf_center,
                                       get_par().epsilon, 
                                       get_par().derotang)); 

    popfirst!(Parameters);
end

function TransRotate(A::AffineTransform2D{Float64}, EPSILON, ANGLE, CENTER, NEWCENTER)
	return translate(rotate(translate(CENTER[1]-EPSILON[1],CENTER[2]-EPSILON[2], A), ANGLE),-NEWCENTER[1], -NEWCENTER[2])
end

function bbox_size(inp_dims::NTuple{2,Integer},
                  A::AffineTransform2D{Float64})
                    
	xmin=typemax(Float64);
	ymin=typemax(Float64);
	xmax=typemin(Float64);
	ymax=typemin(Float64);
	
	xmin_int=typemin(Int64);
	ymin_int=typemin(Int64);
	xmax_int=typemin(Int64);
	ymax_int=typemin(Int64);
	
	width=typemin(Float64);
	height=typemin(Float64);
	
	Ind=[repeat(1:inp_dims[1], inner=inp_dims[2]) repeat(1:inp_dims[2], outer=inp_dims[1])]
	for it=1:(inp_dims[1]*inp_dims[2])
		(x,y)=A(Ind[it,1], Ind[it,2]);
		xmin=min(x, xmin);
		ymin=min(y, ymin);
		xmax=max(x, xmax);
		ymax=max(y, ymax);
	end
	xmin_int=floor(Int64,xmin);
	ymin_int=floor(Int64,ymin);
	xmax_int=ceil(Int64,xmax);
	ymax_int=ceil(Int64,ymax);

	width=xmax_int-xmin_int+1+50;
	height=ymax_int-ymin_int+1+50;
    out_dims=(width, height);

	return(out_dims, (xmin_int, xmax_int, ymin_int, ymax_int))
end

function set_fft_op(PSF::AbstractArray{T,2}, PSFCenter::AbstractArray{T,1}) where {T <: AbstractFloat}
 	MapSize=get_par().cols[1:2]
	MapCenter=floor.(MapSize./2).+1
	MAP=zeros(MapSize)
	ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)

	Id = AffineTransform2D{Float64}()
	Recentering=translate(-(MapCenter[1]-PSFCenter[1]), -(MapCenter[2]-PSFCenter[2]), Id)

	LazyAlgebra.apply!(MAP, ker, Recentering, PSF)
	MAP./=sum(MAP)
	F=FFTOperator(MAP)
	FFT=F\Diag(F*ifftshift(MAP)) .*F;
	
	return FFT
end

function vcreate(::Type{LazyAlgebra.Direct}, A::FieldTransformOperator{T},
                 x::AbstractArray{T, 3}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.cols
    @assert A.cols[3] == 4
    Array{T, 2}(undef, A.rows)
end

function vcreate(::Type{LazyAlgebra.Adjoint}, A::FieldTransformOperator{T},
                x::AbstractArray{T, 2}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.rows
    Array{T,3}(undef, A.cols)
end


function apply!(α::Real,
                ::Type{LazyAlgebra.Direct},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,3},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,2}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.cols
    @assert size(dst) == R.rows
    n = R.rows[2]
    @assert iseven(n)
    z_l_star = R.v_l[1] * src[:,:,1];
    z_r_star = R.v_r[1] * src[:,:,1];
    z_l_disk = R.v_l[1] * src[:,:,2] + R.v_l[2] * src[:,:,3] + R.v_l[3] * src[:,:,4];
    z_r_disk = R.v_r[1] * src[:,:,2] + R.v_r[2] * src[:,:,3] + R.v_r[3] * src[:,:,4];
    dst[:, 1:(n÷2)] = R.H_l_disk * z_l_disk + R.H_l_star * z_l_star;
    dst[:, (n÷2)+1:n] = R.H_r_disk * z_r_disk + R.H_r_star * z_r_star;
    return dst
end

function apply!(α::Real,
                ::Type{LazyAlgebra.Adjoint},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,2},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T, 3}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.rows
    @assert size(dst) == R.cols

    n = R.rows[2]
    @assert iseven(n)
    y_l_star = R.H_l_star' * view(src, :, 1:(n÷2))
    y_l_disk = R.H_l_disk' * view(src, :, 1:(n÷2))
    y_r_star = R.H_r_star' * view(src, :, (n÷2)+1:n)
    y_r_disk = R.H_r_disk' * view(src, :, (n÷2)+1:n)
    I_disk = R.v_l[1] * y_l_disk + R.v_r[1] * y_r_disk
    I_star = R.v_l[1] * y_l_star + R.v_r[1] * y_r_star
    Q = R.v_l[2] * y_l_disk + R.v_r[2] * y_r_disk
    U = R.v_l[3] * y_l_disk + R.v_r[3] * y_r_disk
    dst .= cat(I_star, I_disk, Q, U, dims=3);
    return dst
end

function fg!(x::AbstractArray{T,3}, g::AbstractArray{T,3}) where {T<:AbstractFloat}
    local f::Float64 = 0.0;
    vfill!(g,0.0)
    for k in 1:length(dataset)
        if sum(dataset[k].weights) !=0
            f += fg!(x, g, dataset[k])
        end
    end
    return f
end       

function fg!(x::AbstractArray{T, 3}, g::AbstractArray{T,3}, d::data_table{T}) where {T<:AbstractFloat}
    r = d.H*x - d.data;
    wr = d.weights .* r;
    g .+= d.H'*wr;
    f = vdot(Float64, r,wr);
    return Float64(f)/2
end