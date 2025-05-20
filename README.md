# PADI.jl
**PADI** (Reconstruction of High-contrAst Polarized SOurces and Deconvolution for cIrcumstellar Environments)

## Requirements
- [Julia](https://julialang.org/install/) version 1.11.5 or above
- Python might be required if you want to reproduce the paper's figures. Otherwise not. The python requirements are :
  - MatPlotLib
  - Numpy
  - Skimage
  - astropy
- Other Julia dependencies can be found under the Manifest.toml or Project.toml files
  
Please note that except julia, every requirements will be installed when using the "pkg> precompile" command (Usage section).

## Installation

- Using HTTPS
```sh
git clone https://github.com/SpottedBush/PADI.git
```

- Using SSH:
```sh
git clone git@github.com:SpottedBush/PADI.git
```
N.B.: Using SSH requires you to have a registered SSH key linked to your GitHub account.
Move to the PADI folder then, you need to activate PADI environment using:
```julia
pkg> activate .
pkg> precompile
```

Then, you can check the dependencies with:
```julia
pkg> status
```

## Usage
PADI can be applied using:

```julia
using PADI
x=apply_PADI(x0, A, d, μ)
```
where:

-`x0` is the initialization,

-`A` the convolution by the PSF,

-`d` the dataset uncluding data and weights,

-`μ` a vector of regularization hyperparameters.

## Figures Reproduction

In order to reproduce Padi paper's figure 2, you'll need to execute the files in this exact order:
- test/PADI_data_simulation.jl
- test/generate_double_diff.jl
- test/PADI_reconstruction.jl
- test/PADI_primas/PADI_{regul_type}.jl # The three regularisations presented in the paper

For figure 3 :
- pds70_data/generate_double_diff.jl
- pds70_data/generate_pds70_reconstruction.jl
- pds70_data/{regul_type}_pds70_reconstruction.jl # The three regularisations presented in the paper

Then, for both figures, head to `metrics.ipynb` and execute the cell that contains the function `plot_images_comparison` then execute this function with the right paths.

Each required data for the scripts contained in the test folder can be accessed [here](https://drive.google.com/drive/folders/1RcNd5Qh6XD2Trd07VGPWqlphMU699PFi?usp=drive_link).
Those folders must be extracted in root folder.

## Architecture
<pre>
PADI/
│
├── src/                        # Source files for the PADI method
│   ├── datasimul_tools.jl
│   ├── grad_tools.jl
│   ├── loaders.jl
│   ├── metrics.ipynb           # Jupyter notebook for evaluating metrics and making the figures in the paper
│   ├── padi_methods.jl
│   ├── PADI.jl                 # Main PADI module entry point
│   ├── polarimetric_parameters.jl
│   ├── sure_tools.jl
│   └── utils.jl
│
├── test/                       # Unit tests and simulations
│   ├── PADI_primas/            # PRIMA scripts to find the optimal hyperparameters couple considering a certain ground truth 
│   ├── testsuites/             # Test suites for Grad tools and polarimetric parameters methods
│   ├── generate_double_diff.jl # Requires what PADI_data_simulation.jl produces. Outputs the double difference for the contrast k 
│   ├── PADI_data_simulation.jl # Requires data_for_demo fits file (PSF fits/txt, Iu_star.fits, ddit_simulated_data.fits) produces MASK.fits, DATA.fits, WEIGHT.fits, TRUE_convolved.fits and TRUE.fits for the corresponding contrast_list
│   ├── PADI_reconstruction.jl  # Requires what PADI_data_simulation.jl and generate_double_diff.jl produces and do a single PADI reconstruction (depending on its contrast and regul_type)
│   └── SSIM_diff_from_file.jl  # SSIM comparison on each layer of the input fits files (Ip with Ip etc...)
│
├── LICENCE.md
├── Manifest.toml               # Julia package dependency lock file
├── Project.toml                # Julia project dependencies
├── README.md                   # PADI documentation
│
######## GoogleDrive Part, once downloaded, extract it at the root folder ########
│
├── data_for_demo/              # Contains data used by PADI_data_simulation.jl, needs to be downloaded on the GoogleDrive
│   ├── ddit_simulated_data.fits
│   ├── Iu_star.fits
│   ├── Parameters.txt
│   ├── pds70_angles.txt
│   ├── PSF_centers_Airy.txt
│   └── PSF_parametered_Airy.fits
│  
├── pds70_data/                 # Contains data to reproduce the pds70 figures
│   ├── DATA_processed_coro.fits
│   ├── Ditering.txt
│   ├── generate_double_diff.jl
│   ├── instruments_values_with_crosstalk.txt
│   ├── pds70_angles.txt
│   ├── PSF_centers_Airy.txt
│   ├── PSF_parametered_Airy.fits
│   ├── Results_Separable_DoubleDifference.fits
│   ├── RHAPSODIE.fits
│   ├── WEIGHT_processed_coro.fits
│   ├── Parameters.txt
│   ├── reconstruction_scripts/
│   │   ├── generate_pds70_reconstruction.jl
│   │   ├── rd_pds70_reconstruction.jl
│   │   ├── rj_pds70_reconstruction.jl
│   │   └── rs_pds70_reconstruction.jl
│   └── results/
│       ├── PADI_pds70_mixed_disjoint_{λ_1}_{λ_2}.fits
│       ├── PADI_pds70_mixed_joint_{λ}_{α}.fits
│       └── PADI_pds70_mixed_struct_{λ_1}_{λ_2}_{α}.fits
│  
└── test_results/               # Folder that contains every results produced by the scripts in the "test" folder
    ├── contrast_10e-{value}/   # Results that are produced with 10e-{value} contrast
    │   ├── DATA.fits           # Base DATA with 10e-{value} contrast
    │   ├── MASK.fits           # Pixels MASK
    │   ├── WEIGHT.fits         # Pixels WEIGHT
    │   ├── TRUE.fits           # Ground truth we aim for when reconstructing using PADI
    │   ├── TRUE_convolved.fits # Same as TRUE but convolved
    │   ├── RHAPSODIE.fits      # RHAPSODIE reconstruction for this contrast
    │   ├── Results_Separable_DoubleDifference.fits # DoubleDifference reconstruction for this contrast
    │   └── PADI_method_results/ # PADI reconstructions for this contrast
    │       └── max_iter_{max_iter}/ # Contains results of the PADI methods for max_iter iterations
    │
    └── prima/                  # Results produced by prima scripts
        └── contrast_10e-{value}/
            ├── disjoint_regul/
            │   ├── PADI_{λ_1}_{λ_2}.fits
            │   └── ssim.csv    # Contains SSIM comparison between the different hyperparam couples for the different layers (Iu_disk, Ip_disk, θ.)
            ├── joint_regul/
            │   ├── PADI_{λ}_{α}.fits
            │   └── ssim.csv
            └── struct_regul/
                ├── PADI_{λ_1}_{λ_2}_{α}.fits
                └── ssim.csv
<pre>
