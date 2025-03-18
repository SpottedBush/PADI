#
# sure_tools.jl
#
# Provides an implementation of the Generalized Stein 
# Umbiased Risk Estimator (GSURE) and of the prevision error.
#
# [Stein, 1981] Stein, C. M. (1981). Estimation of the mean of a  
# multivariate normal distribution. The annals of Statistics, 
# pages 1135–1151.
#
# [Eldar, 2008] Eldar, Y. C. (2008). Generalized sure for
# exponential families : Applications to regularization. 
# IEEE Transactions on Signal Processing, 57(2) :471–481.
#
# [Ramani et al., 2008] Ramani, S., Blu, T., and Unser, M. (2008). 
# Monte-Carlo Sure : A Black-Box Optimization of Regularization 
# Parameters for General Denoising Algorithms. IEEE Trans.
# on Image Process., 17(9) :1540–1554.
#
#
#----------------------------------------------------------
#
# This file is part of PADI
#
#
# Copyright (c) 2024-2025 Vincent Tardieux and Laurence Denneulin (see LICENCE.md)
#
#------------------------------------------------


function MSE_data(x_est::Array{T,N}, x_true::Array{T,N}, d::Array{data_table,1}) where {T <: AbstractFloat,N}
    MSE=0.0;
    n=0;
    for data in d
        res = MSE_data(x_est, x_true, data)
        MSE +=res[1];
        n += res[2];
    end
    return MSE,n
end

function MSE_data(x_est::Array{T,N}, x_true::Array{T,N}, d::data_table) where {T <: AbstractFloat,N}
        minimum(d.weights) >=0 || error("invalid weights");    

        res= d.H*(x_est - x_true);
        MSE = vdot( res, d.weights.*res);
        n= count(d.weights .>0);
    return MSE,n
end

function MSE_object(x_est::PolarimetricMap, x_true::PolarimetricMap)
    MSE = zeros(length(fieldnames(PolarimetricMap)) - 1)
    centers=size(x_est)[1:2]./2
    mask=x_true.I_star .> 0
    for ind in CartesianIndices(mask)
        if (ind[1] - centers[1])^2 + (ind[2] - centers[2])^2 < 5^2
            mask[ind]=0.
        end           
    end
    # n_pixels = sum(get_MASK())
    for (i, attr) in enumerate(fieldnames(PolarimetricMap))
        if i == 1 # Skipping field "parameter_type"
            continue
        end
        if i == 11 # Calculating circular MSE for theta field
            # MSE[i - 1] = rad2deg(vnorm2(angle.(exp.(im*2*(x_est.θ - x_true.θ))).*get_MASK()/2)/n_pixels)
            continue
        end
        MSE[i - 1] = vdot((getfield(x_est, attr) - getfield(x_true, attr)).*mask, (getfield(x_est, attr) - getfield(x_true, attr)).*mask)
        # MSE[i - 1] /= vdot(getfield(x_true, attr), getfield(x_true, attr))
    end
    return MSE
end


function absolute_error(x_est::PolarimetricMap, x_true::PolarimetricMap)
    abs_error = zeros(length(fieldnames(PolarimetricMap)) - 1)
    centers=size(x_est)[1:2]./2
    mask=x_true.I_star .> 0
    for ind in CartesianIndices(mask)
        if (ind[1] - centers[1])^2 + (ind[2] - centers[2])^2 < 5^2
            mask[ind]=0.
        end           
    end
    for (i, attr) in enumerate(fieldnames(PolarimetricMap))
        if i == 1 # Skipping field "parameter_type"
            continue
        end
        if i == 11 # Calculating circular MSE for theta field
            # abs_error[i - 1] = rad2deg(vnorm2(angle.(exp.(im*2*(x_est.θ - x_true.θ))).*get_MASK()/2)/n_pixels)
            continue
        end
        x_est_field=getfield(x_est, attr) .* mask
        x_true_field=getfield(x_true, attr) .* mask
        abs_error[i - 1] = sum(abs.(x_est_field[x_true_field .!= 0] - x_true_field[x_true_field .!= 0]))
        #abs_error[i - 1] /= sum(abs.(true[true .!= 0]))
    end
    return abs_error
end

function SSIM(x_est::PolarimetricMap, x_true::PolarimetricMap)
    ssim_values = zeros(length(fieldnames(PolarimetricMap)) - 1)
    for (i, attr) in enumerate(fieldnames(PolarimetricMap))
        if i == 1 # Skipping field "parameter_type"
            continue
        end
        ssim_values[i - 1] = assess_ssim(getfield(x_est, attr), getfield(x_true, attr))
    end
    return ssim_values
end


function MSE_object(x_est::PolarimetricMap, x_true::PolarimetricMap)
    MSE = zeros(length(fieldnames(PolarimetricMap)) - 1)
    n_pixels = 1
    for (i, attr) in enumerate(fieldnames(PolarimetricMap))
        if i == 1 # Skipping field "parameter_type"
            continue
        end
        if i == 7 # Calculating circular MSE for theta field
            MSE[i - 1] = rad2deg(vnorm2(angle.(exp.(im*2*(x_est.θ - x_true.θ))).*1/2)/n_pixels)
            continue
        end
        MSE[i - 1] = vdot(getfield(x_est, attr) - getfield(x_true, attr), getfield(x_est, attr) - getfield(x_true, attr))
        MSE[i - 1] /= vdot(getfield(x_true, attr), getfield(x_true, attr))
    end
    return MSE
end

function sure_crit(x::Array{T,N},
                   δx::Array{T,N}, 
                   d::Array{data_table,1}, 
                   δd::Array{data_table,1}; 
                   mask::Array{Array{Bool,2},1}=Array{Array{Bool,2},1}()) where {T <: AbstractFloat,N}

    data_fidelity=0.0;
    trace= 0.0;
    n=0;

    for k=1:length(d)
        if !isempty(mask)       
            res = sure_crit(x, δx, d[k], δd[k], mask=mask[k])
        else
            res = sure_crit(x, δx, d[k], δd[k])   
        end
        data_fidelity += res[1];
        trace += res[2];
        n += res[3];        
    end
    return [data_fidelity ; trace; n]

end

function sure_crit(x::Array{T,N},
                   δx::Array{T,N}, 
                   d::data_table, 
                   δd::data_table; 
                   mask::Array{Bool,2}=Array{Bool,2}(undef, 0,0)) where {T <: AbstractFloat,N}

        MAP = (d.weights .!=0); #Valid pixels map
        
        !isempty(mask) && (size(MAP)==size(mask)) && (MAP .*= mask);
                
        n=count(MAP);
        if n != 0;
            εb=(δd.data -d.data)[MAP]; #Perturbation of variance ε^2
            εd=d.data - d.H*x; #difference between data and model
            εx=(d.H*(δx - x))[MAP]; #difference between model estimated on perturbed data and on data

            data_fidelity = vdot(εd, d.weights.*εd);  
            
            trace = n*vdot(εb, εx)/vdot(εb, εb)
        else
            data_fidelity=0.
            trace=0.
        end      
    return [data_fidelity ; trace; n]

end

function SURE(solver::Function, 
              par::Array{Float64,1}, 
              d::Array{data_table,1}, 
              δd::Array{data_table,1},       
              x0::Array{Float64,3}; 
              mask::Array{Array{Bool,2},1}=Array{Array{Bool,2},1}())

    x  = solver(x0, d, par);
    δx = solver(x0, δd, par);
 
    return sure_crit(x,δx, d, δd, mask=mask), x
end

function do_data_perturbation(dset::Array{data_table,1},eps::Float64)
    p = data_table[];
    for d in dset
        push!(p,do_data_perturbation(d,eps));
    end
    return p
end                      

function do_data_perturbation(d::data_table, eps::Float64)
    return data_table(d.data + eps* randn(size(d.data)), d.weights, d.H);
end                      


function sure_optim(solver, x)#; linear=true)
        res=sure_tools.SURE(solver, x, grad_tools.dataset, data_perturbed, X0nl)[1]
    return res[1] + 2*res[2];
end
