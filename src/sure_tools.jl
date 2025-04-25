#
# sure_tools.jl
#
# Contains methods to compute indexes.
#
#----------------------------------------------------------
#
# This file is part of PADI
#
#
# Copyright (c) 2024-2025 Vincent Tardieux and Laurence Denneulin (see LICENCE.md)
#
#------------------------------------------------

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