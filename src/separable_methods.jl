#
# separable_methods.jl
#
# Provides the state-of-the-art methods Double Ratio 
# and Double Difference and separable inverse approach.
#
# State of the arts methods have been implemented from:
# [Tinbergen, 2005]  Tinbergen,  J.  (2005).Astronomical Polarimetry. 
# Cambridge  University
#
# Inverse approch :
# [Denneulin, 2020] Approche inverse pour la reconstruction 
# des environnements circumstellaires en polarimétrie 
# avec l’instrument d’imagerie directe ESO/VLT SPHERE IRDIS.. 
# Instrumentation et méthodes pour l'astrophysique [astro-ph.IM]. 
# Université Claude Bernard Lyon 1 (UCBL), 2020. Français. ⟨tel-03200282⟩
#
#
#------------------------------------------------------
#
# This file is part of Rhapsodie
#
# Copyright (c) 2017-2021 Laurence DENNEULIN (see Licence.md)
#TODO : Add a way to return the CRLB


#------------------------------------------------

"""
    Double_Difference(d,ind) -> X
    Double_Difference(d) -> X #if indices loaded with 
    
where X is:
    - a PolarimetricMap if d is of size (N1,N2,K,2) 
    - a PolarimetricPixel if d is of size (K,2).
"""
function Double_Difference(data::Array{Float64,4},ind::Array{Int64,2})
        n1,n2,n3,n4= size(data);
        S=Array{PolarimetricPixel}(undef,n1,n2);
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1   
                S[i1,i2]=Double_Difference(data[i1,i2,:,:], ind);
            end
        end
        return PolarimetricMap(S)
end

function Double_Difference(data::Array{Float64,2},ind::Array{Int64,2})
		Q=mean((data[ind[1,:],1] .-data[ind[1,:],2] 
		          .-(data[ind[2,:],1] .-data[ind[2,:],2]))/2);		
		U=mean((data[ind[3,:],1] .-data[ind[3,:],2] 
		          .-(data[ind[4,:],1] .-data[ind[4,:],2]))/2);
		Iq=(data[ind[1,:],1] .+data[ind[1,:],2] 
		    .+data[ind[2,:],1] .+data[ind[2,:],2])/2;
		Iu=(data[ind[3,:],1] .+data[ind[3,:],2] 
		    .+data[ind[4,:],1] .+data[ind[4,:],2])/2;

        return PolarimetricPixel("stokes", 
                                 (mean(Iq)+ mean(Iu))/2, 
                                 Q,
                                 U)
end

    Double_Difference(data::Array{Float64,4}) = (Double_Difference(data,
                                                                   get_par().indices))

    Double_Difference(data::Array{Float64,2}) = (Double_Difference(data,
                                                                   get_par().indices))

#
#-----------------------------------------------------
#
"""
    Double_Ratio(d,ind) -> X
    Double_Ratio(d) -> X #if indices loaded with 
    
    
where X is:
    - a PolarimetricMap if d is of size (N1,N2,K,2) 
    - a PolarimetricPixel if d is of size (K,2).

"""
function Double_Ratio(data::Array{Float64,4},ind::Array{Int64,2})
        n1,n2,n3,n4= size(data);
        S=Array{PolarimetricPixel}(undef,n1,n2);
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1   
                S[i1,i2]=Double_Ratio(data[i1,i2,:,:], ind);
            end
        end
        return PolarimetricMap(S)
end

function Double_Ratio(data::Array{Float64,2},ind::Array{Int64,2})
		Rq=(data[ind[1,:],1] ./data[ind[1,:],2]) ./(data[ind[2,:],1] ./data[ind[2,:],2]);
		Ru=(data[ind[3,:],1] ./data[ind[3,:],2]) ./(data[ind[4,:],1] ./data[ind[4,:],2]);
		for i=1:length(Rq)
			Rq[i]=sqrt(Rq[i]*(Rq[i]>0));
			Ru[i]=sqrt(Ru[i]*(Ru[i]>0));
		end
		pq=(Rq .-1)./(Rq .+1);
		pu=(Ru .-1)./(Ru .+1);
		Iq=(data[ind[1,:],1] .+data[ind[1,:],2] .+data[ind[2,:],1] .+data[ind[2,:],2])/2;
		Iu=(data[ind[3,:],1] .+data[ind[3,:],2] .+data[ind[4,:],1] .+data[ind[4,:],2])/2;
        
        return PolarimetricPixel("stokes", 
                                 (mean(Iq)+ mean(Iu))/2, 
                                 mean(pq .*Iq),
                                 mean(pu .*Iu))

end


Double_Ratio(data::Array{Float64,4}) = (Double_Ratio(data, get_par().indices))

Double_Ratio(data::Array{Float64,2}) = (Double_Ratio(data, get_par().indices))

#
#-----------------------------------------------------
#
"""
    Linear_Method(d,w) -> X
    Linear_Method(d) -> X #i.e. W=Id
    
    
where X is:
    - a PolarimetricMap if d is of size (N1,N2,K,2) 
    - a PolarimetricPixel if d is of size (K,2).

"""
function Linear_Method(data::Array{Float64,4}, weight::Array{Float64,4})
        n1,n2,n3,n4= size(data);
        x=Array{PolarimetricPixel}(undef,n1,n2);
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1
                if weight[i1,i2] !=0 
                    x[i1,i2]=Linear_Method(data[i1,i2,:,:],weight[i1,i2,:,:]);
                else
                    x[i1,i2] = PolarimetricPixel("stokes",0.0, 0.0, 0.0);
                end
            end
        end
        return PolarimetricMap(x)
end

function Linear_Method(data::Array{Float64,2}, weight::Array{Float64,2})
    T = Float64 # floating point type used for computations

    # Local variables to integrate the normal equations and other
    # quantities.  Only the lower triangular part of the left-hand-side
    # matrix is needed.
    K = get_par().dataset_length # index range for 2nd dimension of `ind`
    @assert size(data)[1] == K;
    @assert size(weight)[1] == K;
    
    local A11::T = 0
    local A21::T = 0
    local A22::T = 0
    local A31::T = 0
    local A32::T = 0
    local A33::T = 0
    local b1::T = 0.0
    local b2::T = 0.0
    local b3::T = 0.0
    local W::T, Wd::T
    for j in 1:2
        W = 0.0
        Wd = 0.0
        for k=1:K
            v1 = get_par().v[k][j][1]
            v2 = get_par().v[k][j][2]
            v3 = get_par().v[k][j][3]
            w, d = weight[k,j], data[k,j]
            W = w
            Wd = w*d

        # Integrate the upper-triangular part of A.
        A11 += v1*W*v1
        A21 += v2*W*v1
        A22 += v2*W*v2
        A31 += v3*W*v1
        A32 += v3*W*v2
        A33 += v3*W*v3

        # Integrate the right-hand side vector b.
        b1 += v1*Wd
        b2 += v2*Wd
        b3 += v3*Wd
        end
    end

    # Solve the normal equations using fast StaticArrays.
    A = SMatrix{3,3,T}((A11, A21, A31,
                        A21, A22, A32,
                        A31, A32, A33))
    b = SVector{3,T}((b1,b2,b3))
    
    #CLRB=inv(A);
    x = A\b
    #x=CLRB*b;
    return PolarimetricPixel("stokes", x)#, CLRB 
end

Linear_Method(data::Array{Float64,4}) = (Linear_Method(data, ones(size(data))))

Linear_Method(data::Array{Float64,2}) = (Linear_Method(data,ones(size(data))))


#
#-----------------------------------------------------
#

"""
    NonLinear_Method(d,w) -> X
    NonLinear_Method(d) -> X #i.e. W=Id
    
    
where X is:
    - a PolarimetricMap if d is of size (N1,N2,K,2) 
    - a PolarimetricPixel if d is of size (K,2).

"""

function NonLinear_Method(data::Array{Float64,4}, weight::Array{Float64,4})
        n1,n2,n3,n4= size(data);
        x=Array{PolarimetricPixel}(undef,n1,n2);
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1   
                if weight[i1,i2] !=0 
                    x[i1,i2]=NonLinear_Method(data[i1,i2,:,:],weight[i1,i2,:,:]);
                 else
                    x[i1,i2] = PolarimetricPixel("intensities",0.0, 0.0, 0.0);
                end                   
            end
        end
        return PolarimetricMap(x)
end

function NonLinear_Method(data::A, weight::A) where {T<:AbstractFloat, A<:AbstractArray{T,2}}
	cost(Θ)=opti_nonlinear(Θ,data,weight)[1]
	(Θ_opt,c)=Bradi.minimize(cost,range(-pi/2,stop=pi/2,length=7))
	f,x=opti_nonlinear(Θ_opt,data,weight)
	return x
end


function opti_nonlinear(Θ::T,data::H, weight::H) where {T<:AbstractFloat, H<:AbstractArray{T,2}}
    K = get_par().dataset_length # index range for 2nd dimension of `ind`
    @assert size(data)[1] == K;
    @assert size(weight)[1] == K;

	local s1 =0.0;
	local s2 =0.0;
	local s3 =0.0;
	local b1 =0.0;
	local b2 =0.0;
	local W::T, Wd::T

	for j=1:2
	    W = 0.0
	    for k=1:K
	        v1 = get_par().v[k][j][1]
            v2 = get_par().v[k][j][2]
            v3 = get_par().v[k][j][3]
            
            V1 = v1;
            V2 = (v1 + v2*cos(2*Θ) + v3*cos(2*Θ));
            
            w, d = weight[k,j], data[k,j]
            W = w
            Wd = w*d    		           		        
	        
	        s1 += V1*W*V1;
            s2 += V1*W*V2;
     		s3 += V2*W*V2;
     		b1 += V1*Wd;
            b2 += V2*Wd;
       end
	end

    A = SMatrix{2,2,T}((s1,s2,s2,s3))
    b = SVector{2,T}((b1,b2))
 
	Iu, Ip = A\b;	
	
	if (Iu>=0) && (Ip>=0)
	    x=PolarimetricPixel("intensities",Iu,Ip,Θ);
		return f_nonlinear(x,data,weight), x
	elseif (Iu<0) || (Ip<0)
		AIu=A*[1;0];
		Iu0=(AIu\b)[];
		if Iu>0
		    x1=PolarimetricPixel("intensities",Iu,0.0,Θ);
		else 
	        x1=PolarimetricPixel("intensities",0.0,0.0,Θ);
	    end
	    
		AIp=A*[0;1];
		Ip=(AIp\b)[];
		if Ip > 0
            x2=PolarimetricPixel("intensities",0.0,Ip,Θ);
	    else 
	        x2=PolarimetricPixel("intensities",0.0,0.0,Θ);
	    end
        
        x3 = PolarimetricPixel("intensities",0.0,0.0,Θ);

		C1=f_nonlinear(x1,data,weight);
		C2=f_nonlinear(x2,data,weight);
		C3=f_nonlinear(x3,data,weight);

        if (C1<=C2)&&(C1<=C3)
            return C1,x1
        elseif (C2<=C1)&&(C2<=C3)
            return C2,x2    
        elseif (C3<=C1)&&(C3<=C2)
            return C3,x3
        else
            error("can't find a good configuration")
	    end
	end
end

function f_nonlinear(x::PolarimetricPixel{T}, data::A, weight::A) where {T<:AbstractFloat, A<:AbstractArray{T,2}}
	local f::T =0.0;
	local W::T, Wd::T
	
	K = get_par().dataset_length # index range
    @assert size(data)[1] == K;
    @assert size(weight)[1] == K;

	for j=1:2
	    W=0.0;
	    Wd =0.0;
	    for k=1:K
            w, d = weight[k,j], data[k,j]
            W = w
            Wd = w*d    		           		        
	        v1 = get_par().v[k][j][1]
            v2 = get_par().v[k][j][2]
            v3 = get_par().v[k][j][3]
            
	        m = v1*x.Iu + (v1 + v2*cos(2*x.θ) + v3*cos(2*x.θ))*x.Ip
	        f += w *(m - d)^2
	    end
	end

	return f
end

NonLinear_Method(data::Array{Float64,4}) = (NonLinear_Method(data, ones(size(data))))

NonLinear_Method(data::Array{Float64,2}) = (NonLinear_Method(data,ones(size(data))))



#=
function FDCR_nonlin(RHO,THETA,W,ind)
    (v,v_norm)=V_calc(alpha, psi)
    C=[cos(THETA); sin(THETA)]
    a11=0;
    a22=0;
    a33=0;
    a12=0;
    a13=0;
    a23=0;
    for i=1:4
        for j=1:2
	        prodsca=abs(dot(v[:,(i-1)*2+j],C))^2
	        dprodsca=cos(THETA)*sin(THETA)*(abs(v[2,(i-1)*2+j])^2-abs(v[1,(i-1)*2+j])^2)+cos(2*THETA)*2*real(v[1,(i-1)*2+j]*conj(v[2,(i-1)*2+j]))
	        a11=a11+sum(W[ind[i,:],j]*v_norm[(i-1)*2+j]^2);
	        a22=a22+sum(W[ind[i,:],j]*prodsca^2);
	        a33=a33+sum(W[ind[i,:],j]*dprodsca^2);		
	        a12=a12+sum(W[ind[i,:],j]*v_norm[(i-1)*2+j]*prodsca);
	        a13=a13+sum(W[ind[i,:],j]*v_norm[(i-1)*2+j]*dprodsca);
	        a23=a23+sum(W[ind[i,:],j]*prodsca*dprodsca);
	    end
    end
    I=[a11 a12 RHO*a13; a12 a22 RHO*a23; RHO*a13 RHO*a23 (RHO^2)*a33];
    if det(I) !=0
        F=inv(I)
        return F
    else
    return zeros(size(I));
    end
end  
=#
