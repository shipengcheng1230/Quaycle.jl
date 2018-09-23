# JuEQ.jl Documentation

## Overview
This is a suite for numerically simulating earthquake sequences in [Julia](https://julialang.org/). The purpose of this package is to provide efficient Julia implementations for simulations in the field of earthquake physics.

Features of this package are listed as below:

- Rate-State Friction Law
- Okada's Dislocation Method
- Boundary Element Method (Quasi-dynamic)


Features to be implemented:
- Viscoelastic Relaxation (priority)
- Fully Elastodynamic Effect
- Off-Fault Materials effect


## Installation
Get the latest version with Julia's package manager:

```julia
] add https://github.com/shipengcheng1230/JuEQ.jl
```

To load the package:

```julia
using JuEQ
```

## Acknowledgements

The simulation of episodic seismic and slow slip events using boundary-element-method is largely benifit from [Yajing Liu](https://liumcgill.wordpress.com/) original Fortran code.
