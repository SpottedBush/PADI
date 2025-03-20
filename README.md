# PADI.jl
**PADI** (Reconstruction of High-contrAst Polarized SOurces and Deconvolution for cIrcumstellar Environments)


## Installation

In the package manager:
- Using HTTPS
```julia
pkg> add https://github.com/SpottedBush/PADI.jl
```

- Using SSH:
```julia
pkg> add git@github.com:SpottedBush/PADI.git
```

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
