module TPolarimetricTestSuite

    using Revise
    using PADI
    using Test

    @testset "TPolarimetricPixel constructor stokes" begin
        pix = TPolarimetricPixel("stokes", 1.0, 2.0, 3.0, 4.0)
        @test pix.I_star == 1.0
        @test pix.I_disk == 2.0
        @test pix.Q == 3.0
        @test pix.U == 4.0
        @test pix.I == pix.I_disk + pix.I_star
        @test pix.I == pix.Iu + pix.Ip_disk
        @test pix.Iu == pix.Iu_disk + pix.Iu_star
        @test pix.Ip_disk == round(sqrt(pix.Q^2 + pix.U^2), digits=3)
        @test pix.I_disk == pix.Iu_disk + pix.Ip_disk
        @test pix.I_star == pix.Iu_star
        @test pix.θ % π == (atan(pix.U, pix.Q) / 2) % π 
    end

    @testset "TPolarimetricPixel constructor intensities" begin
        pix = TPolarimetricPixel("intensities", 1.0, 2.0, 3.0, 4.0)
        @test pix.Iu_star == 1.0
        @test pix.Iu_disk == 2.0
        @test pix.Ip_disk == 3.0
        @test pix.θ == 4.0
        @test pix.I == pix.I_disk + pix.I_star
        @test pix.I == pix.Iu + pix.Ip_disk
        @test pix.Iu == pix.Iu_disk + pix.Iu_star
        @test pix.Ip_disk == round(sqrt(pix.Q^2 + pix.U^2), digits=3)
        @test pix.I_disk == pix.Iu_disk + pix.Ip_disk
        @test pix.I_star == pix.Iu_star
        @test round(pix.θ % π, digits=5) == round((atan(pix.U, pix.Q) / 2) % π, digits=5)
    end

    @testset "TPolarimetricPixel constructor mixed" begin
        pix = TPolarimetricPixel("mixed", 1.0, 2.0, 3.0, 4.0)
        @test pix.Iu_star == 1.0
        @test pix.Iu_disk == 2.0
        @test pix.Q == 3.0
        @test pix.U == 4.0
        @test pix.I == pix.I_disk + pix.I_star
        @test pix.I == pix.Iu + pix.Ip_disk
        @test pix.Iu == pix.Iu_disk + pix.Iu_star
        @test pix.Ip_disk == round(sqrt(pix.Q^2 + pix.U^2), digits=3)
        @test pix.I_disk == pix.Iu_disk + pix.Ip_disk
        @test pix.I_star == pix.Iu_star
        @test pix.θ % π == (atan(pix.U, pix.Q) / 2) % π
    end

    @testset "TPolarimetricPixel operators" begin
        pix_stokes = TPolarimetricPixel("stokes", 1.0, 2.0, 3.0, 4.0)
        pix_intensities = TPolarimetricPixel("intensities", 1.0, -3.0, 5.0, 0.4636476090008061)
        pix_mixed = TPolarimetricPixel("mixed", 1.0, -3.0, 3.0, 4.0)
        @test pix_stokes == pix_intensities
        @test pix_intensities == pix_mixed
        @test pix_stokes == pix_mixed
        @test pix_stokes + pix_stokes == TPolarimetricPixel("stokes", 2.0, 4.0, 6.0, 8.0)
        @test pix_stokes - pix_stokes == TPolarimetricPixel("stokes", 0.0, 0.0, 0.0, 0.0)
    end

    @testset "TPolarimetricMap constructor" begin
        pix = TPolarimetricPixel("stokes", 1.0, 2.0, 3.0, 4.0)
        diff_pix = TPolarimetricPixel("intensities", 1.0, 2.0, 3.0, 4.0)
        pix_list = fill(pix, (2,3))
        
        @test TPolarimetricMap(pix_list) == TPolarimetricMap(pix_list)
        pix_list[1] = diff_pix
        @test_throws "Polarimetric types must match." TPolarimetricMap(pix_list)
    end
end # module