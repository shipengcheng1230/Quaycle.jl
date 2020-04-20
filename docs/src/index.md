# Quaycle.jl Documentation

## Overview
This is a suite for numerically simulating earthquake sequences in [Julia](https://julialang.org/). The purpose of this package is to provide efficient Julia implementations for simulations in the field of earthquake physics.

Features of this package currently:

- ✅ Rate-State Friction
- ✅ Plastic Deformation
- ✅ Boundary Element Method (Quasi-dynamic)
- ✅ Viscoelastic Relaxation
- ✅ Integration with [Gmsh](http://gmsh.info/) and [Paraview](https://www.paraview.org/)

Features to be implemented:
- ❌ Fully elastodynamics
- ❌ Finite element method via [Fenics](https://fenicsproject.org/)

## Installation
Get the latest version with Julia's package manager:

```julia
(v1.3) pkg> add https://github.com/shipengcheng1230/Quaycle.jl
```

To load the package:

```julia
using Quaycle
```
