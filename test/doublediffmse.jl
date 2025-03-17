using PADI
using EasyFITS

const I_HEADER_POS = 1
const I_STAR_HEADER_POS = 2
const I_DISK_HEADER_POS = 3
const Q_HEADER_POS = 4
const U_HEADER_POS = 5
const IU_HEADER_POS = 6
const IP_DISK_HEADER_POS = 9
const THETA_HEADER_POS = 10

true_file = readfits("test_results/contrast_10e-2.0/TRUE.fits")
true_polar_map = PolarimetricMap("mixed",
                                view(true_file, :, :, I_STAR_HEADER_POS)' .+ view(true_file, :, :, I_DISK_HEADER_POS)',
                                view(true_file, :, :, Q_HEADER_POS)',
                                view(true_file, :, :, U_HEADER_POS)', 
                                view(true_file, :, :, IU_HEADER_POS)',
                                view(true_file, :, :, IP_DISK_HEADER_POS)',
                                view(true_file, :, :, THETA_HEADER_POS)')

double_difference_file = readfits("test_results/contrast_10e-2.0/Results_Separable_DoubleDifference.fits")
double_diff = PolarimetricMap("mixed",
            view(double_difference_file, :, :, 1)',
            view(double_difference_file, :, :, 2)',
            view(double_difference_file, :, :, 3)', 
            view(double_difference_file, :, :, 4)',
            view(double_difference_file, :, :, 5)',
            view(double_difference_file, :, :, 6)')

mse = PADI.MSE_object(double_diff, true_polar_map)
println("MSE: ", mse)
