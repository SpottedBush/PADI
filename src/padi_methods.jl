#
# grad_tools.jl
#
# Provide the RHAPSODIE methods to reconstruct the polarimetric 
# parameters from a dataset using:
#
# - VMLM-B (when "mixed" polarimetric parameters are used)
# [Thiébaut, 2002] Thiébaut, E. (2002). Optimization issues in 
# blind deconvolution algorithms. In Astronomical Data Analysis II,
# volume 4847, pages 174–183. International Society for Optics and Photonics.
#
# - Forward-Backward with backtracking (when "stokes" parameters are used)
# [Beck and Teboulle, 2009] Beck, A. and Teboulle, M. (2009). 
# A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems. 
# SIAM J. Imaging Sci., 2(1) :183–202. (TODO: clean implementation)
#
#
# ------------------------------------------------
#
# This file is part of Rhapsodie
#
#
# Copyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

#------------------------------------------------
struct hyperparameters{T<: AbstractFloat}
    λ::T # Weight of the regularization
    ρ::T # Ratio between the weight and the threshold μ = λ/(2*ρ);
end

"""
    apply_PADI('x0', 'A' , 'data', 'hyperparameters'; 'kwd')

where : 
    - x0 is an initial reconstruction
    - A is the operator of the convolution by the PSF
    - par is a set of hyper parameters given as follow : λ_I, λ_Q+U, ρ_I, ρ_Q+U
      where ρ is the ratio between the contribution λ and the threshold μ = λ/2ρ of the hypersmooth regularization;

"""

function apply_PADI(x0::PolarimetricMap, A::D, d::Array{data_table,1}, par::Array{T,1}; mem=3, maxeval=50, maxiter=50, α::Real=1, xtol=(1e-3,1e-8), gtol=(1e-3,1e-8), ftol=(1e-3,1e-8), regul_type::String="struct", verbose::Bool=false) where {T <: AbstractFloat, D <:Mapping}
    n1,n2 = size(x0)
    parameter_type = x0.parameter_type
    X0 = convert(Array{T,3}, x0);
    μ=[hyperparameters(par[1], par[2]); # Iu_star
       hyperparameters(par[3], par[4]); # Iu_disk
       hyperparameters(par[5], par[6])]; # Ip_disk
    lower_born=vcreate(X0);
    upper_born=vcreate(X0);
    if parameter_type == "intensities"
        vfill!(view(lower_born,:,:,1:3),0.0)
        vfill!(view(lower_born,:,:,4),-π)
        vfill!(view(upper_born,:,:,1:3),Inf)
        vfill!(view(upper_born,:,:,4),π)
    elseif parameter_type == "mixed"
        vfill!(view(lower_born,:,:,1:2),0.0)
        vfill!(view(lower_born,:,:,3:4),-Inf)
        vfill!(view(upper_born,:,:,1:4),Inf)
    end
    g=vcreate(X0);
    PADI_fg!(x,g) = apply_gradient!(PolarimetricMap(parameter_type, x), A, g, d, μ, α, regul_type)
    x = vmlmb(PADI_fg!, X0, mem=mem, maxeval=maxeval, maxiter=maxiter, lower=lower_born, upper=upper_born, xtol=xtol,  gtol=gtol, ftol=ftol, verb=verbose);
    return PolarimetricMap(x0.parameter_type, x)
end

"""
    apply_gradient!(X::PolarimetricMap, A::D, g::Array{T,3}, d::Array{data_table,1}, μ::Array{hyperparameters{T},1}, α::Real, regul_type::String) where {T <: AbstractFloat, D <:Mapping}

Applies the gradient to the polarimetric map `X` using the mapping `A`, gradient array `g`, data array `d`, hyperparameters `μ`, regularization parameter `α`, and regularization type `regul_type`.

# Arguments
- `X::PolarimetricMap`: The polarimetric map.
- `A::D`: The mapping.
- `g::Array{T,3}`: The gradient array.
- `d::Array{data_table,1}`: The data array.
- `μ::Array{hyperparameters{T},1}`: The hyperparameters.
- `α::Real`: The regularization parameter.
- `regul_type::String`: The type of regularization ("struct", "joint", or "disjoint").

# Returns
- `f::Float64`: The computed cost function value.
"""
function apply_gradient!(X::PolarimetricMap, A::D, g::Array{T,3}, d::Array{data_table,1}, μ::Array{hyperparameters{T},1}, α::Real, regul_type::String) where {T <: AbstractFloat, D <:Mapping}

    n1, n2, n3 = size(g)
    @assert (n1,n2) == size(X)
    @assert n3 == 4
    
    Ay = cat(X.I_star[:,:], A*X.I_disk[:,:], A*X.Q[:,:], A*X.U[:,:], dims=3)
    # Compute data fidelity term and gradient. (As gradient is initially set to
    # zero, we can recycle it between x and y.)
    @assert size(g) == size(Ay)
    vfill!(g, 0)
    local f::Float64 = 0.0;
    for k = 1:length(d)
        f += fg!(Ay, g, d[k])
    end

    # Convert gradient w.r.t. y into gradient w.r.t. x.  Nothing has to be done
    # for the 2nd and 3rd fields (Q and U) or if Ip = 0.
    @assert size(g) == (n1,n2,4)
    @inbounds for i3 in 2:4
        g[:,:,i3].= A'*view(g,:,:,i3)[:,:];
    end

    if X.parameter_type == "intensities" # Basis under the form (Iu_star, Iu_disk, Ip, θ)
        @inbounds for i2 in 1:n2
            for i1 in 1:n1
                if X.Ip_disk[i1,i2] > 0
                    curr_g_3 = copy(g[i1, i2, 3])
                    g[i1, i2, 3] = g[i1, i2, 2] + cos(2*X.θ[i1, i2]) * g[i1, i2, 3] + sin(2*X.θ[i1, i2]) * g[i1, i2, 4]
                    g[i1, i2, 4] =  -2 * X.Ip_disk[i1, i2] * sin(2*X.θ[i1, i2]) * curr_g_3
                                    + 2 * X.Ip_disk[i1, i2] * cos(2*X.θ[i1, i2]) * g[i1, i2, 4]
                end
            end
        end 
 	    #f+=cost!(μ[1][2] , μ[1][1], X.Iu[:,:], view(g,:,:,1), false);
 	    f+=apply_tikhonov!(X.Iu_star[:,:], view(g,:,:,1), μ[1].λ / (2 * μ[1].ρ));
        f+=apply_edge_preserving_smoothing!(X.Iu_disk[:,:], view(g,:,:,2), μ[2].λ, μ[2].ρ; α)
        f+=apply_edge_preserving_smoothing!(X.Ip_disk[:,:], view(g,:,:,3), μ[3].λ, μ[3].ρ; α)

    elseif X.parameter_type == "mixed" # Basis under the form (Iu_star, Iu_disk, Q, U)
        @inbounds for i2 in 1:n2
            for i1 in 1:n1
                if X.Ip_disk[i1,i2] > 0
                    g[i1, i2, 3] += g[i1, i2, 2] * (X.Q[i1, i2] / X.Ip_disk[i1, i2])
                    g[i1, i2, 4] += g[i1, i2, 2] * (X.U[i1, i2] / X.Ip_disk[i1, i2])
                end
            end
        end
        
 	    #f+=cost!(μ[1][2] , μ[1][1], X.Iu[:,:], view(g,:,:,1), false);

        # Struct regularization
        if regul_type == "struct"
            tmp_grad = zeros(T, n1, n2, 2)
            f+=apply_tikhonov!(X.Iu_star[:,:], view(g,:,:,1), μ[1].λ / (2 * μ[1].ρ));
            f+=apply_edge_preserving_smoothing!(cat(X.Iu_disk[:,:], X.Ip_disk[:,:], dims=3), tmp_grad, μ[2].λ, μ[2].ρ; α=α^2)
            f+=apply_struct_regul!(X.Iu_disk, view(tmp_grad,:,:,1), μ[3].λ * α)
            f+=apply_struct_regul!(X.Ip_disk, view(tmp_grad,:,:,2), μ[3].λ)
            g[:,:,2] .+= tmp_grad[:,:,1]
            g[:,:,3] .+= X.Q .* tmp_grad[:,:,2] ./ X.Ip_disk
            g[:,:,4] .+= X.U .* tmp_grad[:,:,2] ./ X.Ip_disk

        # Joint regularization
        elseif regul_type == "joint"
            tmp_grad = zeros(T, n1, n2, 2)
            f+=apply_tikhonov!(X.Iu_star[:,:], view(g,:,:,1), μ[1].λ / (2 * μ[1].ρ));
            f+=apply_edge_preserving_smoothing!(cat(X.Iu_disk[:,:], X.Ip_disk[:,:], dims=3), tmp_grad, μ[2].λ, μ[2].ρ; α=α)
            g[:,:,2] .+= tmp_grad[:,:,1]
            g[:,:,3] .+= X.Q .* tmp_grad[:,:,2] ./ X.Ip_disk
            g[:,:,4] .+= X.U .* tmp_grad[:,:,2] ./ X.Ip_disk
        
        # Disjoint regularization
        elseif regul_type == "disjoint"
            tmp_grad = zeros(T, n1, n2)
            f+=apply_tikhonov!(X.Iu_star[:,:], view(g,:,:,1), μ[1].λ / (2 * μ[1].ρ));
            f+=apply_edge_preserving_smoothing!(X.Iu_disk[:,:], view(g,:,:,2), μ[2].λ, μ[2].ρ)
            f+=apply_edge_preserving_smoothing!(X.Ip_disk[:,:], tmp_grad, μ[3].λ, μ[3].ρ)
            g[:,:,3] .+= X.Q .* tmp_grad ./ X.Ip_disk
            g[:,:,4] .+= X.U .* tmp_grad ./ X.Ip_disk
        end

    elseif X.parameter_type == "stokes" # Basis under the form (Iu_star, Iu_disk, Q, U)
 	    f+=apply_tikhonov!(X.I_star[:,:], view(g,:,:,1), μ[1].λ / (2 * μ[1].ρ));
        f+=apply_edge_preserving_smoothing!(cat(X.I_disk[:,:], X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:4), μ[2].λ, μ[2].ρ; α)
 	    #f+=cost!(μ[1][2] , μ[1][1], X.I[:,:], view(g,:,:,1), false);    
     end
	#f+=cost!(μ[2][2] , μ[2][1], cat(X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:3), false);

	return f
end

"""
    apply_struct_regul!(x::AbstractArray{T,2}, g::AbstractArray{T,2}, λ::Real) where {T <: AbstractFloat}

Applies structural regularization to the array `x` and updates the gradient array `g`.

# Arguments
- `x::AbstractArray{T,2}`: The input array.
- `g::AbstractArray{T,2}`: The gradient array.
- `λ::Real`: The regularization parameter.

# Returns
- `sum(abs, x) * λ`: The computed cost function value.
"""
function apply_struct_regul!(x::AbstractArray{T,2},
        g::AbstractArray{T,2},
        λ::Real) where {T <: AbstractFloat}
    for i in 1:size(x, 1)
        for j in 1:size(x, 2)
            g[i, j] += λ
        end
    end
    return sum(abs, x) * λ
end

"""
    apply_tikhonov!(x::AbstractArray{T,2}, g::AbstractArray{T,2}, λ::Real) where {T <: AbstractFloat}

Applies Tikhonov regularization to the array `x` and updates the gradient array `g`.

# Arguments
- `x::AbstractArray{T,2}`: The input array.
- `g::AbstractArray{T,2}`: The gradient array.
- `λ::Real`: The regularization parameter.

# Returns
- `f::T`: The computed cost function value.
"""
function apply_tikhonov!(x::AbstractArray{T,2},
        g::AbstractArray{T,2},
        λ::Real) where {T <: AbstractFloat}
    m,n = size(x)                               
    f = zero(T);
    r = zero(T);
    x1 = zero(T);
    x2 = zero(T);

    for i=1:m-1
        for j=1:n-1
            x1= (x[i,j] - x[i+1,j])/2
            x2= (x[i,j] - x[i,j+1])/2
            r = x1^2 + x2^2;
            ## Cost functon ##
            f += λ * r / 2;
                ## Gradient in x ##
                g[i,j] += λ*(x1 + x2);
                g[i+1,j] -= λ*x1; 
                g[i,j+1] -= λ*x2; 
        end
    end
    return f
end

"""
    apply_edge_preserving_smoothing!(x::AbstractArray{T,3}, g::AbstractArray{T,3}, λ::Real, ρ::Real; α=1.0) where {T <: AbstractFloat}

Applies edge-preserving smoothing to the 3D array `x` and updates the gradient array `g`.

# Arguments
- `x::AbstractArray{T,3}`: The input array.
- `g::AbstractArray{T,3}`: The gradient array.
- `λ::Real`: The regularization parameter.
- `ρ::Real`: The edge-preserving parameter.
- `α::Real`: The regularization parameter for the first dimension (default is 1.0).

# Returns
- `f::T`: The computed cost function value.
"""
function apply_edge_preserving_smoothing!(x::AbstractArray{T,3},
    g::AbstractArray{T,3},
    λ::Real,
    ρ::Real;
    α=1.0) where {T <: AbstractFloat}

    m,n,o = size(x)
    f = zero(T);
    r = zero(T);
    μ = λ/(2*ρ);
    x1 = zero(T);
    x2 = zero(T);
        
    for i=1:m-1
        for j=1:n-1
            ndx=zero(T)
            for k=1:o
                x1= (x[i,j,k] - x[i+1,j,k])/2
                x2= (x[i,j,k] - x[i,j+1,k])/2
                
                nrm = x1^2 + x2^2;
                (k==1) && (nrm *= α) 
                ndx += nrm
            end
            r =  ndx + μ^2;
            ## Cost functon ##
            f += λ*(√r -  μ);
            if r>0
                for k=1:o
                    x1= (x[i,j,k] - x[i+1,j,k])/2
                    x2= (x[i,j,k] - x[i,j+1,k])/2

                    (k==1) && (x1 *= α)
                    (k==1) && (x2 *= α)
                    ## Gradient in x ##
                    ∂r=2*√r;
                    g[i,j,k] += λ*(x1 + x2)/∂r;
                    g[i+1,j,k] -= λ*x1/∂r; 
                    g[i,j+1,k] -= λ*x2/∂r; 
                    
                end
            end
        end
    end
    
    return f
end

"""
    apply_edge_preserving_smoothing!(x::AbstractArray{T,2}, g::AbstractArray{T,2}, λ::Real, ρ::Real) where {T <: AbstractFloat}

Applies edge-preserving smoothing to the 2D array `x` and updates the gradient array `g`.

# Arguments
- `x::AbstractArray{T,2}`: The input array.
- `g::AbstractArray{T,2}`: The gradient array.
- `λ::Real`: The regularization parameter.
- `ρ::Real`: The edge-preserving parameter.

# Returns
- `f::T`: The computed cost function value.
"""
function apply_edge_preserving_smoothing!(x::AbstractArray{T,2},
    g::AbstractArray{T,2},
    λ::Real,
    ρ::Real) where {T <: AbstractFloat}

    m, n = size(x)
    f = zero(T)
    r = zero(T)
    μ = λ / (2 * ρ)
    x1 = zero(T)
    x2 = zero(T)

    for i = 1:m-1
        for j = 1:n-1
            x1 = (x[i, j] - x[i+1, j]) / 2
            x2 = (x[i, j] - x[i, j+1]) / 2

            nrm = x1^2 + x2^2

            r = nrm + μ^2
            ## Cost function ##
            f += λ * (√r - μ)
            if r > 0
                ## Gradient in x ##
                ∂r = 2 * √r
                g[i, j] += λ * (x1 + x2) / ∂r
                g[i+1, j] -= λ * x1 / ∂r
                g[i, j+1] -= λ * x2 / ∂r
            end
        end
    end
    return f
end
