# PADI.jl
**PADI** (Reconstruction of High-contrAst Polarized SOurces and Deconvolution for cIrcumstellar Environments)


## Installation

- Using HTTPS
```sh
git clone https://github.com/SpottedBush/PADI.jl
```

- Using SSH:
```sh
git clone git@github.com:SpottedBush/PADI.git
```
N.B.: Using SSH requires you to have a registered SSH key linked to your GitHub account.

## Usage
Move to the PADI folder then, you need to activate PADI environment using:
```julia
pkg> activate .
pkg> precompile
```

Then, you can check the dependencies with:
```julia
pkg>status
```
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


Further usage can be found in the test folder where reconstruction scripts examples are given.
Each required data for the scripts contained in the test folder can be accessed here : https://drive.google.com/drive/folders/1RcNd5Qh6XD2Trd07VGPWqlphMU699PFi?usp=drive_link

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
│   ├── PADI_data_simulation.jl # Requires data_for_demo fits file (PSF fits/txt, Iu_star.fits, ddit_simulated_data.fits) produces MASK.fits, DATA.fits, WEIGHT.fits, TRUE_convolved.fits and TRUE.fits for the corresponding contrast_list
│   ├── PADI_reconstruction.jl  # Requires what PADI_data_simulation.jl produces and do a single PADI reconstruction (depending on its contrast and regul_type)
│   └── SSIM_diff_from_file.jl  # SSIM comparison on each layer of the input fits files (Ip with Ip etc...).
│
├── LICENCE.md
├── Manifest.toml               # Julia package dependency lock file
├── Project.toml                # Julia project dependencies
└── README.md                   # PADI documentation
<pre>
