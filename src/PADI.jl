#
# PADI.jl
#
# Package for the Polarimetric and Angular Differential Imaging of high-contrast circumstellar environments (PADI)
#
#----------------------------------------------------------
#
# This file is part of PADI
#
#
# Copyright (c) 2024-2025 Vincent Tardieux and Laurence Denneulin (see LICENCE.md)
#
#------------------------------------------------

module PADI

    export
        PolarimetricPixel,
        PolarimetricMap,
        write_polar_map,
        read_and_fill_polar_map,
        convert,
        load_parameters,
        get_par,
        load_data,
        fg!,
        SetCropOperator,
        crop,
        crop!,
        pad,
        set_fft_op,
        TwoDimensionalTransformInterpolator,
        FieldTransformOperator,
        data_simulator,
        generate_model,
        data_generator,
        generate_parameters,
        apply_PADI,
        apply_edge_preserving_smoothing!,
        SSIM

    import Base: +, -, *, /, ==, getindex, setindex!, read, write, convert

    using OptimPackNextGen
    import OptimPackNextGen: BraDi #va devenir BraDi avec un D majuscule
    using SpecialFunctions
    using TwoDimensional
    using FFTW
    using LinearInterpolators
    using Statistics
    using LinearAlgebra
    using LazyAlgebra
    import LazyAlgebra: Mapping, vcreate, vcopy, apply!
    using StaticArrays
    using FITSIO
    using EasyFITS
    using DelimitedFiles
    using Random
    using ImageQualityIndexes
    import ImageQualityIndexes: assess_ssim


    include("polarimetric_parameters.jl")
    include("datasimul_tools.jl")
    include("grad_tools.jl")
    include("loaders.jl")
    include("padi_methods.jl")
    include("sure_tools.jl")
    include("utils.jl")
end