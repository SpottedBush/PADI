module TGradToolsTestSuite
    using Revise
    using PADI
    using Test
    using LazyAlgebra
    using TwoDimensional
    using LinearInterpolators
    
    @testset "FieldTransformOperator_simple_case" begin
        T = FieldTransformOperator((2, 2, 4), (2, 4), (1.0, 0.0, 1.0), (0.0, 1.0, 0.0), LazyAlgebra.Id, LazyAlgebra.Id, LazyAlgebra.Id, LazyAlgebra.Id)
        for i in range(1, 10)
            x = randn(Float64, 2, 2, 4)
            y = randn(Float64, 2, 4)
            @test vdot(T * x, y) / vdot(x, T' * y) ≈ 1.0
        end
    end
    @testset "FieldTransformOperator_translating" begin
	    Id = TwoDimensional.AffineTransform2D{Float64}()
        ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)
        translate1 = TwoDimensional.translate(3e-1, 8e-1, Id)
        translate2 = TwoDimensional.translate(0.2, 1.4, Id)
        input_size = (2, 2)
        output_size = (2, 2)
        T_l_star = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, translate1)
        T_l_disk = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, translate2)
        T_r_star = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, translate1)
        T_r_disk = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, translate2)
        
        T = FieldTransformOperator((2, 2, 4), (2, 4), (1.0, 0.0, 1.0), (0.0, 1.0, 0.0), T_l_star, T_l_disk, T_r_star, T_r_disk)
        for i in range(1, 10)
            x = randn(Float64, 2, 2, 4)
            y = randn(Float64, 2, 4)
            @test vdot(T * x, y) / vdot(x, T' * y) ≈ 1.0
        end
    end
end