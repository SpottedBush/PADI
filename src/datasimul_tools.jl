#
# datasimul_tools.jl
#
# Provide tools to simulate synthetic parameters and dataset.
#
#----------------------------------------------------------
#
# This file is part of PADI
#
#
# Copyright (c) 2024-2025 Vincent Tardieux and Laurence Denneulin (see LICENCE.md)
#
#------------------------------------------------

"""
	data_generator(model::AbstractArray{T,N}, weights::AbstractArray{T,N}; bad=zero(T)) where {T<:AbstractFloat, N}

Generate data based on the given model and weights. The function adds noise to the model data based on the weights.

# Arguments
- `model::AbstractArray{T,N}`: The model data array.
- `weights::AbstractArray{T,N}`: The weights array.
- `bad`: The value to assign when the weight is zero (default is `zero(T)`).

# Returns
- `data`: The generated data array with added noise.

# Throws
- `error("invalid weights")` if any weight is not finite or is negative.
"""
function data_generator(model::AbstractArray{T,N}, weights::AbstractArray{T,N};bad=zero(T)) where {T<:AbstractFloat,N}   
    #seed === nothing ||  Random.seed!(seed);
    
    data = Array{T}(undef, size(model));
    @inbounds for i in eachindex(data, weights)
        w = weights[i]
        (isfinite(w) && w >= 0 ) || error("invalid weights")
        if w > 0            
            data[i] = model[i] + randn()/sqrt(w)    
        elseif w == 0 
            data[i] = bad;
        end
    end
    return data
end

"""
	generate_model(S::PolarimetricMap, A::Mapping)

Generate a model based on the given polarimetric map and mapping.

# Arguments
- `S::PolarimetricMap`: The input polarimetric map.
- `A::Mapping`: The mapping to apply.

# Returns
- `M`: The generated model array.
"""
function generate_model(S::PolarimetricMap, A::Mapping)
    @assert size(S) == get_par().cols[1:2];
    
    M=Array{Float64,3}(undef, get_par().rows[1], 
                              get_par().rows[2], 
                              get_par().dataset_length)
	ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)    
    for k=1:get_par().dataset_length
        output_size = (get_par().rows[1], Int64(get_par().rows[2]/2));
        input_size = get_par().cols[1:2];
    	T_star_left = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][1])
    	T_disk_left = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][2])
    	T_star_right = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][3])
    	T_disk_right = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][4])
	    F = FieldTransformOperator(get_par().cols, 
	                                        get_par().rows, 
	                                        get_par().v[k][1],
	                                        get_par().v[k][2],
	                                        T_star_left,
											T_disk_left,
	                                        T_star_right,
											T_disk_right)
	    
	    M[:,:,k].= F * cat(S.I_star, A * S.I_disk, A * S.Q, A * S.U, dims=3);
	end
    return M
end
    
"""
	tdata_simulator(Good_Pix, tau, A::Mapping; ro_noise=8.5)

Simulate data based on the given parameters.

# Arguments
- `Good_Pix`: The array indicating good pixels.
- `tau`: The polarisation ratio of the disk. (tau = ip/(iu+ip)
- `A::Mapping`: The mapping to apply.
- `ro_noise`: The noise level (default is 8.5).

# Returns
- `D`: The simulated data array.
- `W`: The weights array.
- `S`: The generated polarimetric map.
- `CS`: The convolved polarimetric map.
"""
function tdata_simulator(Good_Pix, tau, A::Mapping; ro_noise=8.5)
    map_size=get_par().cols[1:2];
    
    S = generate_parameters(map_size, tau);

    M = generate_model(S, A);
    
    VAR = max.(M,zero(eltype(M))) .+ro_noise^2
	W = Good_Pix ./ VAR
	D = data_generator(M, W)
	
	check_MSE(M, D, W);
	
    CS = PolarimetricMap("stokes", S.I_star, A*S.I_disk, A*S.Q, A*S.U)
	return D, W, S, CS
end

"""
	ddit_data_simulator(Good_Pix, A::Mapping, S::PolarimetricMap; ro_noise=8.5)

Simulates data for a given polarimetric map.

# Arguments
- `Good_Pix`: A matrix indicating good pixels.
- `A::Mapping`: A mapping object.
- `S::PolarimetricMap`: A polarimetric map object.
- `ro_noise`: (Optional) The noise level, default is 8.5.

# Returns
- `D`: The generated data.
- `W`: The weight matrix.
- `S`: The padded polarimetric map.
- `CS`: The corrected polarimetric map.

# Notes
- If the size of the polarimetric map `S` is not equal to the expected size, it will be padded.
- The function generates a model `M` based on the polarimetric map `S` and mapping `A`.
- The variance `VAR` is calculated based on the model `M` and noise level `ro_noise`.
- The function adjusts `Good_Pix` based on the centers and a threshold.
- The weight matrix `W` is calculated as the ratio of `Good_Pix` to `VAR`.
- The data `D` is generated using the model `M` and weight matrix `W`.
- The function checks the mean squared error (MSE) of the model `M`, data `D`, and weight matrix `W`.
- The function returns the generated data `D`, weight matrix `W`, padded polarimetric map `S`, and corrected polarimetric map `CS`.
"""
function ddit_data_simulator(Good_Pix, A::Mapping, S::PolarimetricMap; ro_noise=8.5)
    if size(S) != get_par().cols[1:2]
        @warn "Size of the Polarimetric Map is wrong, it will be padded."
        S = pad(S);
    end
   
    M = generate_model(S, A);
    
    VAR = max.(M, zero(eltype(M))) .+ ro_noise^2
	centers=[get_par().center, get_par().center - [get_par().epsilon[1][2][1], get_par().epsilon[1][2][2]] + [0, get_par().rows[2]/2]]
    for ind in CartesianIndices(Good_Pix)
        if (ind[1] - centers[1][1])^2 + (ind[2] - centers[1][2])^2 < 15^2
            Good_Pix[ind]=0.
        end
		if (ind[1] - centers[2][1])^2 + (ind[2] - centers[2][2])^2 < 15^2
            Good_Pix[ind]=0.
        end
    end
	W = Good_Pix ./ VAR
	D = data_generator(M, W)
	check_MSE(M, D, W);
    S = PolarimetricMap("stokes", S.I_star, S.I_disk, S.Q, S.U)
	CS = PolarimetricMap("stokes", S.I_star, A*S.I_disk, A*S.Q, A*S.U)
	return D, W, S, CS
end

"""
	check_MSE(model, data, weights)

Check the mean squared error (MSE) between the model and data, weighted by the weights.

# Arguments
- `model`: The model data array.
- `data`: The generated data array.
- `weights`: The weights array.

# Prints
- The MSE, the number of valid weights, and the MSE per valid weight.
"""
function check_MSE(model, data, weights)
	MSE = vdot(data - model, weights .* (data-model)) ;
	N=count(weights .> 0);
	println("MSE=$MSE, N=$N, MSE/N=$(MSE/N)");
end

"""
	generate_parameters(map_size, tau)

Generate the parameters for the polarimetric map based on the given map size and polarisation ratio of the disk.

# Arguments
- `map_size`: The size of the map.
- `tau`: The polarisation ratio of the disk. (tau = ip/(iu+ip)

# Returns
- `PolarimetricMap`: The generated polarimetric map with intensities.
"""
function generate_parameters(map_size, tau)
	Ip=zeros(map_size);
	Iu=zeros(map_size);
	θ=zeros(map_size);
	STAR1=zeros(map_size);
	STAR2=zeros(map_size);
	
	for i=1:map_size[1]
    	for j=1:map_size[2]
    		r1=(map_size[1]+1)/2-i;
    		r2=(map_size[2]+1)/2-j;
    		if (r1^2+r2^2<=20^2)
        		Iu[i,j]=1000;
        		Ip[i,j]=tau*Iu[i,j]/(1-tau);
    		end 
    		if ((r1^2+r2^2>=25^2)&&(r1^2+r2^2<=27^2))
        		Iu[i,j]=1000;
        		Ip[i,j]=tau*Iu[i,j]/(1-tau);
    		end
    		if ((r1^2+r2^2>=32^2)&&(r1^2+r2^2<=40^2))
        		Iu[i,j]=1000;
        		Ip[i,j]=tau*Iu[i,j]/(1-tau);
    		end
    		θ[i,map_size[2]+1-j]=atan(j-map_size[2]/2,i-map_size[1]/2);
			STAR1[i,j]=200*exp(-((i-map_size[1]/2)^2+(j-map_size[2]/2)^2)/(2*75^2))
			STAR2[i,j]=100000*exp(-((i-map_size[1]/2)^2+(j-map_size[2]/2)^2)/(2*7^2))
			if ((((map_size[1]+1)/2-i)^2+((map_size[2]+1)/2-j)^2)<=10^2)
        		STAR2[i,j]=800;
        		Iu[i,j]=0;
        		Ip[i,j]=0;	
    		end
			if ((((map_size[1]+1)/2-i)^2+((map_size[2]+1)/2-j)^2)<=70^2)
        		STAR1[i,j]=50;		
    		end
		end
	end    
	θ = θ.*(Ip.!=0);
	STAR = STAR1+STAR2
	STAR[round(Int64,10*map_size[1]/16)-3,round(Int64,10*map_size[2]/16)] = 20000.0;
	STAR[round(Int64,10*map_size[1]/16),round(Int64,10*map_size[2]/16)-3] = 100000.0;

    return PolarimetricMap("intensities", STAR, Iu, Ip, θ);
end
