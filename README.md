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

First, you need to activate PADI environment using:
```julia
pkg> activate .julia/packages/PADI/XXXXX
pkg> precompile
```
where XXXXX is the folder version.

You can check the dependencies with:

```julia
pkg>status
```
PADI can be applied using:

```julia
x=apply_PADI(x0, A, d, μ)
```
where:

-`x0` is the initialization,

-`A` the convolution by the PSF,

-`d` the dataset uncluding data and weights,

-`μ` a vector of regularization hyperparameters.


Further usage can be found in the test folder where reconstruction scripts examples are given.
Each required data for the scripts contained in the test folder can be accessed here : https://drive.google.com/drive/folders/1RcNd5Qh6XD2Trd07VGPWqlphMU699PFi?usp=drive_link
