function SetCropOperator(table::Vector{NTuple{N, AffineTransform2D}}) where {N}
    DSIZE=get_par().rows[1]
    MASK=ones(get_par().cols[1:2]);
    for k=1:length(table)
    X1=table[k][1](1,1);
    X2=table[k][1](1,DSIZE);
    X3=table[k][1](DSIZE, DSIZE);
    X4=table[k][1](DSIZE, 1);
    Mask1=zeros(get_par().cols[1:2]);

    for i=1:get_par().cols[1]
	    for j=1:get_par().cols[2]	
		    if ((i-X1[1])*(X2[1]-X1[1]) +(j-X1[2])*(X2[2]-X1[2]) >0)&&((i-X2[1])*(X3[1]-X2[1]) +(j-X2[2])*(X3[2]-X2[2]) >0)&&((i-X3[1])*(X4[1]-X3[1]) +(j-X3[2])*(X4[2]-X3[2]) >0) &&((i-X4[1])*(X1[1]-X4[1]) +(j-X4[2])*(X1[2]-X4[2]) >0)	
		    Mask1[i,j,:] .=1;
		    end
	    end
    end
    X1=table[k][2](1,1);
    X2=table[k][2](1,DSIZE);
    X3=table[k][2](DSIZE, DSIZE);
    X4=table[k][2](DSIZE, 1);
    Mask2=zeros(get_par().cols[1:2]);

    for i=1:get_par().cols[1]
	    for j=1:get_par().cols[2]	
		    if ((i-X1[1])*(X2[1]-X1[1]) +(j-X1[2])*(X2[2]-X1[2]) >0)&&((i-X2[1])*(X3[1]-X2[1]) +(j-X2[2])*(X3[2]-X2[2]) >0)&&((i-X3[1])*(X4[1]-X3[1]) +(j-X3[2])*(X4[2]-X3[2]) >0) &&((i-X4[1])*(X1[1]-X4[1]) +(j-X4[2])*(X1[2]-X4[2]) >0)	
		    Mask2[i,j,:] .=1;
		    end
	    end
    end
    MASK .*= Mask1.*Mask2
    end
    #write(FitsFile, "MASK.fits", MASK, overwrite=true)
    push!(MASK_save, MASK)
end

function crop(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    #@assert size(X) .==   get_par().cols
    xymin_1=findfirst(get_MASK() .!= 0.)
    xymin_2=findfirst(get_MASK()' .!= 0.)
    xymax_1=findlast(get_MASK() .!= 0.)
    xymax_2=findlast(get_MASK()' .!= 0.)
    
    xmin=min(xymin_1[1], xymin_2[2]);
    ymin=min(xymin_1[2], xymin_2[1]);
    xmax=max(xymax_1[1], xymax_2[2]);
    ymax=max(xymax_1[2], xymax_2[1]);
    X .*= get_MASK()
    Y=X[xmin:xmax, ymin:ymax];
    Y[.!isfinite.(Y)].=0
    return Y
end

function crop!(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    #@assert size(X) .==   get_par().cols 
    X[.!isfinite.(X)].=0
    X .*= get_MASK()
end
        
function crop(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,3}}
    #@assert size(X) .==   get_par().cols
    Y=copy(X);
    return crop!(Y)
end

function crop!(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,3}}
    #@assert size(X) .==   get_par().cols
    for k in size(X)[3] 
        crop!(view(X,:,:,k))
    end     
end

function crop(X::TPolarimetricMap{T}) where {T<:AbstractFloat}
    #@assert size(X) .==   get_par().cols
    return TPolarimetricMap(X.parameter_type,
                           crop!(view(X.I,:,:)),
                           crop!(view(X.I_star,:,:)),
                           crop!(view(X.I_disk,:,:)),
                           crop!(view(X.Q,:,:)),
                           crop!(view(X.U,:,:)),
                           crop!(view(X.Iu,:,:)),
                           crop!(view(X.Iu_star,:,:)),
                           crop!(view(X.Iu_disk,:,:)),
                           crop!(view(X.Ip_disk,:,:)),        
                           crop!(view(X.θ,:,:)))   
end        

function crop!(X::TPolarimetricMap{T})  where {T<:AbstractFloat}
        crop!(view(X.I,:,:));
        crop!(view(X.I_star,:,:));
        crop!(view(X.I_disk,:,:));
        crop!(view(X.Q,:,:));
        crop!(view(X.U,:,:));
        crop!(view(X.Iu,:,:));
        crop!(view(X.Iu_star,:,:));
        crop!(view(X.Iu_disk,:,:));
        crop!(view(X.Ip_disk,:,:));        
        crop!(view(X.θ,:,:));
end        
        
function pad(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    X_size = size(X);
    center_diff = X_size./2 .- get_par().cols[1:2]./2;
    ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)    
    Id = AffineTransform2D{Float64}()
    center_change = translate(center_diff[1],center_diff[2], Id)   
    PAD=TwoDimensionalTransformInterpolator(get_par().cols[1:2], X_size, ker, ker, center_change)
    return PAD*X;
end         
        
function pad(X::TPolarimetricMap{T}) where {T<:AbstractFloat}  
    #@assert size(X) .==   get_par().cols
    return TPolarimetricMap(X.parameter_type,
                           pad(view(X.I,:,:)),
                           pad(view(X.I_star,:,:)),
                           pad(view(X.I_disk,:,:)),
                           pad(view(X.Q,:,:)),
                           pad(view(X.U,:,:)),
                           pad(view(X.Iu,:,:)),
                           pad(view(X.Iu_star,:,:)),
                           pad(view(X.Iu_disk,:,:)),
                           pad(view(X.Ip_disk,:,:)),        
                           pad(view(X.θ,:,:)))   

end