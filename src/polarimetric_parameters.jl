#
# polarimetric_parameters.jl
#
# Contains the pixel and map structs. 
#
#----------------------------------------------------------
#
# This file is part of PADI
#
#
# Copyright (c) 2024-2025 Vincent Tardieux and Laurence Denneulin (see LICENCE.md)
#
#------------------------------------------------

using EasyFITS
#------------------------------------------------
# PADI Struct definition
# For pixels
struct PolarimetricPixel{T<: AbstractFloat} 
    parameter_type::AbstractString   #either "stokes", "intensities" or "mixed"
    I::T
    I_star::T
    I_disk::T
    Q::T
    U::T
    Iu::T
    Iu_star::T
    Iu_disk::T
    Ip_disk::T # We supposed that the star's light is unpolarized
    θ::T
end

#For map of pixels
struct PolarimetricMap{T<: AbstractFloat} 
    parameter_type::AbstractString   #either "stokes", "intensities" or "mixed"
    I::Array{T, 2};
    I_star::Array{T, 2};
    I_disk::Array{T, 2};
    Q::Array{T, 2};
    U::Array{T, 2};
    Iu::Array{T, 2};
    Iu_star::Array{T, 2};
    Iu_disk::Array{T, 2};
    Ip_disk::Array{T, 2}; # We supposed that the star's light is unpolarized
    θ::Array{T, 2};
end


PolarimetricMap(parameter_type::AbstractString,        
                    I::AbstractArray{T, 2},
                    I_star::AbstractArray{T, 2},
                    I_disk::AbstractArray{T, 2},
                    Q::AbstractArray{T, 2},
                    U::AbstractArray{T, 2},
                    Iu::AbstractArray{T, 2},
                    Iu_star::AbstractArray{T, 2},
                    Iu_disk::AbstractArray{T, 2},
                    Ip_disk::AbstractArray{T, 2},
                    θ::AbstractArray{T, 2}) where {T<:AbstractFloat} = 
                        PolarimetricMap(parameter_type,        
                        convert(Array{T},I),
                        convert(Array{T},I_star),
                        convert(Array{T},I_disk),
                        convert(Array{T},Q),
                        convert(Array{T},U),
                        convert(Array{T},Iu),
                        convert(Array{T},Iu_star),
                        convert(Array{T},Iu_disk),
                        convert(Array{T},Ip_disk),
                        convert(Array{T},θ))
#------------------------------------------------
# Constructors 
"""
    PolarimetricPixel(parameter_type, x) -> PolarimetricPixel
    PolarimetricMap(parameter_type, x) -> PolarimetricMap
    
create an object of type PolarimetricParameter from either:
    - Parameters I_star, I_disk, Q, U (i.e. parameter_type = 'stokes')
    - Parameters Iu_star, Iu_disk, Ip and θ (i.e. parameter_type = 'intensities')
    - Parameters Iu_star, Iu_disk, Q, U (i.e. parameter_type = 'mixed')
    
Each parameter can be called from the structure. For exemple with a 
construction from Stokes parameters S=(I_star, I_disk, Q, U):

    using PADI
    X = PolarimetricParameter(S, 'stokes');
    X.I_disk #yields the Stokes parameter I_disk
    X.Ip_disk #yields the polarized intensity Ip_disk
    X[1,1] #yields a PolarimetricPix at the CartesianIndex (1,1);
    X[1,1].I_disk #yields the Stokes parameter I_disk at the CartesianIndex (1,1); 

    PolarimetricMap(parameter_type, n1, n2) -> PolarimetricMap
    
yields an empty    
"""    
    function PolarimetricPixel(parameter_type::AbstractString, 
                             x1::T,
                             x2::T, 
                             x3::T,
                             x4::T) where {T<:AbstractFloat}
        if parameter_type == "stokes"
            I_star = x1            # Stokes parameter I_star
            I_disk = x2            # Stokes parameter I_disk
            Q = x3            # Stokes parameter Q
            U = x4            # Stokes parameter U
            Ip_disk = sqrt(Q^2 + U^2)    # intensity of linearly polarized light
            Iu = I_star + I_disk - Ip_disk
            Iu_star = I_star 
            Iu_disk = I_disk - Ip_disk
            θ = atan(U, Q) / 2 # angle of linearly polarized light
            I = I_star + I_disk
        elseif parameter_type == "intensities"
            Iu_star = x1
            I_star = x1           # Hypothesis: Star's light is unpolarized
            Iu_disk = x2
            Ip_disk = x3
            I_disk = Iu_disk + Ip_disk
            Iu = Iu_star + Iu_disk
            θ = x4       # angle of linearly polarized light
            I = Iu_star + Iu_disk + Ip_disk   # Stokes parameter I (total light intensity)
            Q = Ip_disk * cos(2*θ);    # Stokes parameter Q
            U = Ip_disk * sin(2*θ);    # Stokes parameter U
        elseif parameter_type == "mixed"
            Iu_star = x1
            I_star = Iu_star
            Iu_disk = x2
            Iu = Iu_star + Iu_disk
            Q = x3            # Stokes parameter Q
            U = x4            # Stokes parameter U  
            Ip_disk = sqrt(Q^2 + U^2)    # intensity of linearly polarized light
            I_disk = Iu_disk + Ip_disk
            I = Iu_star + Iu_disk + Ip_disk     # Stokes parameter I (total light intensity)
            θ = atan(U, Q)/2 # angle of linearly polarized light
        else
            error("unkown type, known types are : 'stokes', 'intensities' and 'mixed'")
        end
        PolarimetricPixel(parameter_type, I, I_star, I_disk, Q, U, Iu, Iu_star, Iu_disk, Ip_disk, θ)
    end

    function PolarimetricPixel(parameter_type::AbstractString, 
                             x::N) where {N<:AbstractArray{Float64, 1}}
        @assert length(x) == 4
        PolarimetricPixel(parameter_type,
                        x[1],
                        x[2],
                        x[3],
                        x[4]);
    end 

    function PolarimetricMap(x::Matrix{PolarimetricPixel{T}}) where {T <: AbstractFloat}
        n1, n2 = size(x)
        par_type=x[1].parameter_type;          
        I = Array{T}(undef, n1, n2)    # Stokes parameter I (total light intensity)
        I_star = Array{T}(undef, n1, n2)
        I_disk = Array{T}(undef, n1, n2)
        Q = Array{T}(undef, n1, n2)    # Stokes parameter Q
        U = Array{T}(undef, n1, n2)    # Stokes parameter U
        Iu = Array{T}(undef, n1, n2)
        Iu_star = Array{T}(undef, n1, n2)
        Iu_disk = Array{T}(undef, n1, n2)
        Ip_disk = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
        θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1   
                pix= x[i1,i2];
                if par_type == pix.parameter_type
                   I[i1, i2] = pix.I;
                   I_star[i1, i2] = pix.I_star;
                   I_disk[i1, i2] = pix.I_disk;
                   Q[i1, i2] = pix.Q;
                   U[i1, i2] = pix.U;
                   Iu_star[i1, i2] = pix.Iu_star;
                   Iu_disk[i1, i2] = pix.Iu_disk;
                   Iu[i1, i2] = pix.Iu;
                   Ip_disk[i1, i2] = pix.Ip_disk;
                   θ[i1, i2] = pix.θ;
                else
                   error("Polarimetric types must match.")
                end
            end
        end
        PolarimetricMap(par_type, I, I_star, I_disk, Q, U, Iu, Iu_star, Iu_disk, Ip_disk, θ)
    end
 
    function PolarimetricMap(parameter_type::AbstractString, 
                             x1::A,
                             x2::A,
                             x3::A,
                             x4::A) where {T<:AbstractFloat, A<:AbstractArray{T,2}}
        n1, n2 = size(x1)
        @assert ((n1,n2) == size(x2)) && ((n1,n2) == size(x3) && ((n1,n2) == size(x4)))
        pixel_list = Matrix{PolarimetricPixel{T}}(undef, n1, n2)
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1
            pix = PolarimetricPixel(parameter_type, x1[i1, i2], x2[i1, i2], x3[i1, i2], x4[i1, i2])
            pixel_list[i1, i2] = pix
            end
        end
        PolarimetricMap(pixel_list)
    end

    function PolarimetricMap(parameter_type::AbstractString, 
                             x::Array{T,3}) where {T<:AbstractFloat}
        n1, n2, n3 = size(x)

        @assert n3 == 4
        PolarimetricMap(parameter_type, 
                        copy(x[:,:,1]),
                        copy(x[:,:,2]),
                        copy(x[:,:,3]),
                        copy(x[:,:,4]));
    end
    
    function PolarimetricMap(parameter_type::AbstractString, n1::Int, n2::Int)
        return PolarimetricMap(parameter_type,
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2))
    end
    
#------------------------------------------------
# Base fonction redefinitions
    
    Base.size(A::PolarimetricMap{T}) where {T<:AbstractFloat} = size(A.I)
    
    Base.length(A::PolarimetricMap{T}) where {T<:AbstractFloat} = prod(size(A)) * 4
    
    Base.length(A::PolarimetricPixel{T}) where {T<:AbstractFloat} = 4
    
    Base.getindex(X::PolarimetricMap{T}, i::CartesianIndex{2}) where {T<:AbstractFloat} =
    PolarimetricPixel(X.parameter_type, X.I[i], X.I_star[i], X.I_disk[i], X.Q[i], X.U[i], X.Iu[i], X.Iu_star[i], X.Iu_disk[i], X.Ip_disk[i], X.θ[i])
    
    Base.getindex(X::PolarimetricMap{T}, i::Int) where {T<:AbstractFloat} = 
    PolarimetricPixel(X.parameter_type, X.I[i], X.I_star[i], X.I_disk[i], X.Q[i], X.U[i], X.Iu[i], X.Iu_star[i], X.Iu_disk[i], X.Ip_disk[i], X.θ[i])
    
    Base.getindex(X::PolarimetricMap{T}, i::Int, j::Int) where {T<:AbstractFloat} = 
    getindex(X, CartesianIndex(i,j))

    function Base.setindex!(X::PolarimetricMap{T}, x::PolarimetricPixel{T}, i::Int64, j::Int64) where {T<:AbstractFloat}
        X.I[i,j]=x.I;
        X.I_star[i,j]=x.I_star;
        X.I_disk[i,j]=x.I_disk;
        X.Q[i,j]=x.Q;
        X.U[i,j]=x.U;
        X.Iu[i,j]=x.Iu;
        X.Iu_star[i,j]=x.Iu_star;
        X.Iu_disk[i,j]=x.Iu_disk;
        X.Ip_disk[i,j]=x.Ip_disk;
        X.θ[i,j]=x.θ;
    end

    function ==(a::PolarimetricMap{T}, b::PolarimetricMap{T}) where {T<:AbstractFloat}
        if a.parameter_type != b.parameter_type
            @warn "Parameter type differs."
        end
        return round.(a.I) == round.(b.I) &&
        round.(a.I_star) == round.(b.I_star) &&
        round.(a.I_disk) == round.(b.I_disk) &&
        round.(a.Q) == round.(b.Q) &&
        round.(a.U) == round.(b.U) &&
        round.(a.Iu) == round.(b.Iu) &&
        round.(a.Iu_star) == round.(b.Iu_star) &&
        round.(a.Iu_disk) == round.(b.Iu_disk) &&
        round.(a.Ip_disk) == round.(b.Ip_disk) &&
        round.(a.θ) == round.(b.θ)
    end

    function ==(a::PolarimetricPixel{T}, b::PolarimetricPixel{T}) where {T<:AbstractFloat}
        if a.parameter_type != b.parameter_type
            @warn "Parameter type differs."
        end
        return round(a.I, digits=5) == round(b.I, digits=5) &&
        round(a.I_star, digits=5) == round(b.I_star, digits=5) &&
        round(a.I_disk, digits=5) == round(b.I_disk, digits=5) &&
        round(a.Q, digits=5) == round(b.Q, digits=5) &&
        round(a.U, digits=5) == round(b.U, digits=5) &&
        round(a.Iu, digits=5) == round(b.Iu, digits=5) &&
        round(a.Iu_star, digits=5) == round(b.Iu_star, digits=5) &&
        round(a.Iu_disk, digits=5) == round(b.Iu_disk, digits=5) &&
        round(a.Ip_disk, digits=5) == round(b.Ip_disk, digits=5) &&
        round(a.θ, digits=5) == round(b.θ, digits=5)
    end

    function +(x::PolarimetricMap{T}, y::PolarimetricMap{T}) where {T<:AbstractFloat}
        I_star = x.I_star + y.I_star
        I_disk = x.I_disk + y.I_disk
        Q = x.Q + y.Q
        U = x.U + y.U
        return PolarimetricMap("stokes", I_star, I_disk, Q, U)
    end
    function +(x::PolarimetricPixel{T}, y::PolarimetricPixel{T}) where {T<:AbstractFloat}
        I_star = x.I_star + y.I_star
        I_disk = x.I_disk + y.I_disk
        Q = x.Q + y.Q
        U = x.U + y.U
        return PolarimetricPixel("stokes", I_star, I_disk, Q, U)
    end
        
    function -(x::PolarimetricPixel{T}, y::PolarimetricPixel{T}) where {T<:AbstractFloat}
        I_star = x.I_star - y.I_star
        I_disk = x.I_disk - y.I_disk
        Q = x.Q - y.Q
        U = x.U - y.U
        return PolarimetricPixel("stokes", I_star, I_disk, Q, U)
    end

    function -(x::PolarimetricMap{T}, y::PolarimetricMap{T}) where {T<:AbstractFloat}
        I_star = x.I_star - y.I_star
        I_disk = x.I_disk - y.I_disk
        Q = x.Q - y.Q
        U = x.U - y.U
        return PolarimetricMap("stokes", I_star, I_disk, Q, U)
    end
      
    vcopy(x::PolarimetricMap{T}) where {T<:AbstractFloat} = PolarimetricMap(x.parameter_type,
                                                 x.I,
                                                 x.I_star,
                                                 x.I_disk,
                                                 x.Q,
                                                 x.U,
                                                 x.Iu,
                                                 x.Iu_star,
                                                 x.Iu_disk,
                                                 x.Ip_disk,
                                                 x.θ)


    function vcreate(x::PolarimetricMap{T}) where {T<:AbstractFloat}
        @assert (x.parameter_type == "stokes") |
                (x.parameter_type == "intensities") |
                (x.parameter_type == "mixed")
        n1,n2=size(x);
        return PolarimetricMap(x.parameter_type, n1, n2)
    end

    function +(x::PolarimetricMap{T}, y::Array{T,3}) where {T<:AbstractFloat} 
        @assert size(y)[1:2] == size(x)       
        if x.parameter_type == "stokes"
           I_star = x.I_star + view(y, :, :, 1);
           I_disk = x.I_disk + view(y, :, :, 2);
           Q = x.Q + view(y, :, :, 3);
           U = x.U + view(y,:,:,4);
           return PolarimetricMap("stokes", I_star, I_disk, Q, U)
        elseif x.parameter_type == "intensities"
           Iu_star = x.Iu_star + view(y, :, :, 1);
           Iu_disk = x.Iu_disk + view(y, :, :, 2);
           Ip_disk = x.Ip_disk + view(y, :, :, 3);
           θ = x.θ + view(y, :, :, 4);
           return PolarimetricMap("intensities", Iu_star, Iu_disk, Ip_disk, θ)
        elseif x.parameter_type == "mixed"
           Iu_star = x.Iu_star + view(y, :, :, 1);
           Iu_disk = x.Iu_disk + view(y, :, :, 2);
           Q = x.Q + view(y, :, :, 3);
           U = x.U + view(y, :, :, 4);
           return PolarimetricMap("mixed", Iu_star, Iu_disk, Q, U)
        else
            error("unknown parameter type")
        end
    end
    
    +(y::Array{T,3}, x::PolarimetricMap{T}) where {T<:AbstractFloat} = x + y
    -(x::PolarimetricMap{T}, y::Array{T,3}) where {T<:AbstractFloat} = x + (-y)
     
     
    function convert(::Type{Array{T,3}}, x::PolarimetricMap{T}) where {T <:AbstractFloat}
         if x.parameter_type == "stokes"
           return cat(x.I_star, x.I_disk, x.Q, x.U, dims=3)
        elseif x.parameter_type == "intensities"
           return cat(x.Iu_star, x.Iu_disk, x.Ip_disk, x.θ, dims=3)
        elseif x.parameter_type == "mixed"
           return cat(x.Iu_star, x.Iu_disk, x.Q, x.U, dims=3)
        else
            error("unknown parameter type")
        end
    end

const I_HEADER_POS = 1
const I_STAR_HEADER_POS = 2
const I_DISK_HEADER_POS = 3
const Q_HEADER_POS = 4
const U_HEADER_POS = 5
const IU_HEADER_POS = 6
const IU_STAR_HEADER_POS = 7
const IU_DISK_HEADER_POS = 8
const IP_DISK_HEADER_POS = 9
const THETA_HEADER_POS = 10
#------------------------------------------------
# Writting function to save PolarimetricMap in fits file
"""
    write(X,'filename.fits') 
    
where X is a PolarimetricMap, write a fitsfile

"""
function write_polar_map(X::PolarimetricMap{T}, filename::AbstractString; overwrite::Bool = false) where {T<:AbstractFloat}
    data = cat(X.I', X.I_star', X.I_disk', X.Q', X.U', X.Iu', X.Iu_star', X.Iu_disk', X.Ip_disk', X.θ', dims=3)
    writefits(filename,
    ["D" => ("Ok", "")],
    data, overwrite=overwrite)
end

"""
    read('parameter_type','filename.fits') -> PolarimetricMap
    
create an object of type PolarimetricMap from a fits file with:
    - Parameters I, Q, U (i.e. parameter_type = 'stokes')
    - Parameters Iu, Ip and θ (i.e. parameter_type = 'intensities')
    - Parameters Iu, Q, U (i.e. parameter_type = 'mixed')
   
"""

function read_and_fill_polar_map(parameter_type::AbstractString, filename::AbstractString)
    @assert (parameter_type == "stokes") |
                (parameter_type == "intensities") |
                (parameter_type == "mixed")
    X=readfits(filename);
    if parameter_type == "stokes"
        return PolarimetricMap(parameter_type,
                                view(X, :, :, I_STAR_HEADER_POS)',
                                view(X, :, :, I_DISK_HEADER_POS)',
                                view(X, :, :, Q_HEADER_POS)',
                                view(X, :, :, U_HEADER_POS)')

    elseif parameter_type == "intensities"
        return PolarimetricMap(parameter_type,
                                view(X, :, :, IU_STAR_HEADER_POS)',
                                view(X, :, :, IU_DISK_HEADER_POS)',
                                view(X, :, :, IP_DISK_HEADER_POS)',
                                view(X, :, :, THETA_HEADER_POS)')

    elseif parameter_type == "mixed"
        return PolarimetricMap(parameter_type,
                                view(X, :, :, IU_STAR_HEADER_POS)',
                                view(X, :, :, IU_DISK_HEADER_POS)',
                                view(X, :, :, Q_HEADER_POS)',
                                view(X, :, :, U_HEADER_POS)')

    else
        return PolarimetricMap(parameter_type,
                           view(X, :, :, I_HEADER_POS)',
                           view(X, :, :, I_STAR_HEADER_POS)',
                           view(X, :, :, I_DISK_HEADER_POS)',
                           view(X, :, :, Q_HEADER_POS)',
                           view(X, :, :, U_HEADER_POS)',
                           view(X, :, :, IU_HEADER_POS)',
                           view(X, :, :, IU_STAR_HEADER_POS)',
                           view(X, :, :, IU_DISK_HEADER_POS)',
                           view(X, :, :, IP_DISK_HEADER_POS)',
                           view(X, :, :, THETA_HEADER_POS)')
    end
end