# Visualize

This package provides some utilities for writing results to VTK (vts, vti, vtu, vtm, pvd, etc...)
  file for postprocessing with [Paraview](https://www.paraview.org/). Data cache for
  writing results are also implemented where users could explore other meshing tools.

Notice that these function are not fully tested. If you encounter any problem,
  please file an issue with your MWE.

## Public Interface
```@autodocs
Modules = [JuEQ]
Pages = ["visualize.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```
