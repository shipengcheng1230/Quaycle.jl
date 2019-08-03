# Visualize

This package provides some utilities for writing results to VTK (vts, vti, vtu, vtm, pvd, etc...)
  file for postprocessing with [Paraview](https://www.paraview.org/). Data cache for
  writing results are also implemented where users could explore other meshing tools.

One short note is that if you store 6 components vector, for example ``σ``,
  [Paraview](https://www.paraview.org/) will interpret them in the order of
  ``σ_{xx}``, ``σ_{yy}``, ``σ_{zz}``, ``σ_{xy}``, ``σ_{yz}``, ``σ_{xz}``.
  Whereas in our [`assemble`](@ref) for stress, the order, same as e.g OpenFOAM, is
  ``σ_{xx}``, ``σ_{xy}``, ``σ_{xz}``, ``σ_{yy}``, ``σ_{yz}``, ``σ_{zz}``.

Notice that these function are not fully tested. If you encounter any problem,
  please file an issue with your MWE.

## Public Interface
```@autodocs
Modules = [Quaycle]
Pages = ["visualize.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```
