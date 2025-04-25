using PADI
using DelimitedFiles
using EasyFITS

# contrast_list = [i for i in range(-1.5, 0, step=0.5)]
contrast_list = [-2.0]
max_iter = 700
# α=10^-5
λ = -2.998282567943823
α= 5.209978160792237

par=readdlm("data_for_demo/Parameters.txt")
DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
# DerotAng = deg2rad.(readdlm("data_for_demo/pds70_angles.txt", Float64)[1:64])
max_angle = 64
DerotAng = [deg2rad(i) for i in range(1, max_angle, length=64)]
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

parameter_type = "mixed"

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

PADI.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)

for k in contrast_list
    println("------Contrast parameter: ", k, "------")
    true_polar_map = PADI.read_and_fill_polar_map(parameter_type, "test_results/contrast_10e$(k)/TRUE.fits")
    diff_fits = EasyFITS.readfits("test_results/contrast_10e$(k)/Results_Separable_DoubleDifference.fits")
    diff_polar_map = PADI.PolarimetricMap(parameter_type, diff_fits[ :, :, 1]', diff_fits[ :, :, 5]', diff_fits[:, :, 2]', diff_fits[:, :, 3]')

    PADI.load_data("test_results/contrast_10e$(k)/DATA.fits", "test_results/contrast_10e$(k)/WEIGHT.fits")

    PSF = readfits("data_for_demo/PSF_parametered_Airy.fits");
    A = set_fft_op(PSF[1:end÷2,:]'[:,:], psf_center[1:2]);
    X0 = diff_polar_map;
    regularisation_parameters = 10 .^[-8, -1., λ, -0.66, -3, -0.66] #(in log10) star, disk | Iu_star, Iu_star, Iu_disk, Iu_disk, Ip_disk, Ip_disk

    x = apply_PADI(X0, A, PADI.dataset, regularisation_parameters, α=10^α, regul_type="joint", maxeval=1000, maxiter=max_iter, verbose=true);
    crop!(x)
    write_polar_map(x, "test_results/contrast_10e$(k)/PADI_method_results/max_iter_$(max_iter)/PADI_opti_params_$(parameter_type)_joint_regul_from_dd_on_$(regularisation_parameters[5]).fits", overwrite=true)
    empty!(PADI.dataset)
end
