# PADI.jl
**PADI** (Reconstruction of High-contrAst Polarized SOurces and Deconvolution for cIrcumstellar Environments)


## Installation

In the package manager:
#TOCHANGE
```julia
pkg> add https://github.com/LaurenceDenneulin/Rhapsodie.jl
```
if you use HTTPS or:

```julia
pkg> add git@github.com:LaurenceDenneulin/Rhapsodie.jl
```
if you use SSH.

## Usage

First, you need to activate Rhapsodie environment using:
```julia
pkg> activate .julia/packages/Rhapsodie/XXXXX
pkg> precompile
```
where XXXXX is the folder version.

You can check the dependencies with:

```julia
pkg>status
```
Rhapsodie can be applied using:

```julia
x=apply_rhapsodie(x0, A, d, μ)
```

where:

-`x0` is the initialization,

-`A` the convolution by the PSF,

-`d` the dataset uncluding data and weights,

-`μ` a vector of regularization hyperparameters.
