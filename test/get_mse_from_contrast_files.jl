using PADI
using DelimitedFiles
using EasyFITS

max_iter = 700
α=10^-5
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

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

PADI.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)
mse_list = Vector{Vector{Float64}}()

for k in range(-3, 0, step=0.5)
    parameter_type = "stokes"
    base_path = joinpath("test_results", "contrast_10e$(k)")

    x_true = read_and_fill_polar_map(parameter_type, joinpath(pwd(), base_path, "TRUE.fits"))
    files = readdir(joinpath(base_path, "PADI_method_results", "max_iter_700"))
    for file in files
        if endswith(file, ".fits")
            parts = split(basename(file), '_')
            lambda = parse(Float64, parts[2])
            alpha = parse(Float64, parts[3][1:end-5])  # remove ".fits" extension

            x_est = read_and_fill_polar_map(parameter_type, joinpath(base_path, "PADI_method_results", "max_iter_700", file))
            res = MSE_object(x_est, x_true)

            Iu_disk_mse = res[8]
            Ip_disk_mse = res[9]
            theta_mse = res[10]

            println("File: ", file, " | Lambda: ", lambda, " | Alpha: ", alpha, " | Iu_disk_mse: ", Iu_disk_mse, " | Ip_disk_mse: ", Ip_disk_mse, " | Theta_mse: ", theta_mse)

            push!(mse_list, [lambda, alpha, k, Iu_disk_mse, Ip_disk_mse, theta_mse])
        end
    end
end

writedlm("test_results/mse_list_test.txt", mse_list)