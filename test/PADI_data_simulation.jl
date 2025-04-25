using PADI
using DelimitedFiles
using EasyFITS
using Dates

if prod(readdir() .!= "test_results")     
    mkdir("test_results")
end

contrast_list = [-2.0]
# contrast_list = [i for i in range(-3, 0, step=0.5)] # decoment this line to test multiple contrasts

DSIZE=256;
NTOT=64;
Nframe=2;
Nrot=NTOT
Nangle=NTOT ÷ (Nframe*4)
Center=[DSIZE/2, DSIZE/2];
# DerotAng = deg2rad.(readdlm("data_for_demo/pds70_angles.txt", Float64)[1:64])
max_angle = 64
DerotAng = [deg2rad(i) for i in range(1, max_angle, length=64)]

Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],[10.7365 , -1.39344]));
end

psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");
PADI.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)

writedlm("data_for_demo/Parameters.txt", [DSIZE; NTOT; Nframe; Nrot; Center; 300; [10.7365 , -1.39344]]);

psf = readfits("data_for_demo/PSF_parametered_Airy.fits");
const A=set_fft_op((psf[1:end÷2,:]'), get_par().psf_center[1]);

Iu_star_fits = readfits("data_for_demo/Iu_star.fits");
Iu_star_original = view(Iu_star_fits, :, :, 1) * 50

ddit_fits = readfits("data_for_demo/ddit_simulated_data.fits");
I_disk_original = view(ddit_fits, :, :, 1)

for k in contrast_list
    Iu_star = copy(Iu_star_original)
    I_disk = copy(I_disk_original)

    contrast = 10.0^k
    normalization_term = contrast * maximum(Iu_star) / maximum(I_disk)

    I_disk .*= normalization_term

    Ip_disk = ddit_fits[:,:,2]
    Ip_disk .*= normalization_term
    scattering = ddit_fits[:,:,3]

    Iu_disk = I_disk - Ip_disk

    Iu_star = Matrix(Iu_star)

    X0 = PADI.PolarimetricMap("intensities", Iu_star, Iu_disk, Ip_disk, scattering)

    GoodPixMap = rand(0.0:1e-16:1.0,(DSIZE, 2*DSIZE)).< 0.99;

    data, weight, S, S_convolved = PADI.ddit_data_simulator(GoodPixMap, A, X0, ro_noise=8.5);

    writefits("test_results/contrast_10e$(k)/MASK.fits", ["D" => ("Ok", "")], PADI.get_MASK()', overwrite=true)

    writefits("test_results/contrast_10e$(k)/DATA.fits",
    ["DATE" => (now(), "date of creation")],
    mapslices(transpose,data,dims=[1,2]), overwrite=true)

    writefits("test_results/contrast_10e$(k)/WEIGHT.fits",
    ["DATE" => (now(), "date of creation")],
    mapslices(transpose,weight,dims=[1,2]), overwrite=true)

    S_convolved = PADI.crop(S_convolved)
    S = PADI.crop(S)

    PADI.write_polar_map(S_convolved, "test_results/contrast_10e$(k)/TRUE_convolved.fits", overwrite=true)
    PADI.write_polar_map(S, "test_results/contrast_10e$(k)/TRUE.fits", overwrite=true)
end