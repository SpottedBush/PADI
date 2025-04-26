using PADI
using PRIMA
using DelimitedFiles
using EasyFITS
import Base.Filesystem: mkpath


par=readdlm("data_for_demo/Parameters.txt")
DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
max_angle = 64
DerotAng = [deg2rad(i) for i in range(1, max_angle, length=64)]
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");
PADI.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)


PSF = readfits("data_for_demo/PSF_parametered_Airy.fits");
A = set_fft_op(PSF[1:end÷2,:]'[:,:],psf_center[1:2]);
header = Vector{Vector{String}}()
push!(header, ["contrast", "λ", "α", "Iu_disk_ssim", "Ip_disk_ssim", "theta_ssim"])

# for k in range(-3, 0, step=0.5)
regul_type = "joint"
k = -2.0
ssim_list = Vector{Vector{Any}}()
root_path = "test_results/contrast_10e$(k)/"
dir_path = "test_results/prima/contrast_10e$(k)"
if !isdir(dir_path)
    mkdir(dir_path)
end
PADI.load_data("$(root_path)DATA.fits", "$(root_path)WEIGHT.fits")

function calculate_ssim_for_prima(X::Vector{Float64})
    λ, α = X
    regularisation_parameters = 10 .^[-8, -1., λ, λ + 3, 0, -1.]
    regularisation_parameters[5] = 0
    println("Iteration : λ: ", λ, " α: ", α)
    diff_fits = EasyFITS.readfits("test_results/contrast_10e$(k)/Results_Separable_DoubleDifference.fits")
    diff_polar_map = PADI.PolarimetricMap("mixed", diff_fits[ :, :, 1]', diff_fits[ :, :, 5]', diff_fits[:, :, 2]', diff_fits[:, :, 3]')
    X0 = diff_polar_map
    x_est = apply_PADI(X0, A, PADI.dataset, regularisation_parameters, maxeval=500, maxiter=1000, α=10^α, regul_type=regul_type)
    true_polar_map = PADI.read_and_fill_polar_map("mixed", "$(root_path)TRUE.fits")
    crop!(x_est)
    mkpath("test_results/prima/contrast_10e$(k)/$(regul_type)_regul/")
    write_polar_map(x_est, "test_results/prima/contrast_10e$(k)/$(regul_type)_regul/PADI_$(λ)_$(α).fits", overwrite=true)
    curr_ssim = PADI.SSIM(x_est, true_polar_map)
    ssim_entry = [k, λ, α, curr_ssim[8], curr_ssim[9], curr_ssim[10]]
    println("SSIM values: Iu_disk: ", curr_ssim[8], " | Iu_star ", curr_ssim[9])
    push!(ssim_list, ssim_entry)
    return 1 - (sum(curr_ssim[8:9]) / 2)
end
optimal_hyperparams, info = PRIMA.newuoa(calculate_ssim_for_prima, [-1., -1.], rhobeg=4, rhoend=1e-1, maxfun=20)
io = open("test_results/prima/contrast_10e$(k)/$(regul_type)_regul/ssim.csv", "w")
writedlm(io, header, ',')
writedlm(io, ssim_list, ',')
close(io)
println("Optimal hyperparameters lambda, alpha: ", optimal_hyperparams[1], " | ", optimal_hyperparams[2])
println("Info: ", info)