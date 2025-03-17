using PADI
using EasyFITS
using DelimitedFiles

parameter_type = "mixed"
prima_path = "test_results/prima/contrast_10e-2.0/" * "joint_regul/"
base_path = "test_results/contrast_10e-2.0/"
true_polar_map = PADI.read_and_fill_polar_map(parameter_type, base_path * "TRUE.fits")

error_list = Vector{Vector{Float64}}()

files = readdir(prima_path)
for file in files
    if endswith(file, ".fits")
        parts = split(basename(file), '_')
        lambda = parse(Float64, parts[2])
        alpha = parse(Float64, parts[3][1:end-5])  # remove ".fits" extension
        
        x_est = read_and_fill_polar_map(parameter_type, prima_path * file)
        res = PADI.SSIM(x_est, true_polar_map)
        # res = PADI.MSE_object(x_est, true_polar_map)

        Iu_disk_error = res[8]
        Ip_disk_error = res[9]
        theta_error = res[10]
        println("File: ", file, " | Lambda: ", lambda, " | Alpha: ", alpha, " | Iu_disk_error: ", Iu_disk_error, " | Ip_disk_error: ", Ip_disk_error, " | Theta_error: ", theta_error)
        push!(error_list, [lambda, alpha, Iu_disk_error, Ip_disk_error, theta_error])
    end
end
writedlm(prima_path * "not_norm_abs_error.csv", error_list, ',')
