# HDF5 Utilities

To use [HDF5](https://github.com/JuliaIO/HDF5.jl) functionality,
```julia
pkg> add HDF5

julia> using Quaycle
julia> using HDF5
```

This package provides a bunch of utilities, such as storing simulation properties
  and writing ODE solution on the fly to HDF5 format file for postprocessing on
  different platform and package version. We currently don't support automatically saving
  green's function to HDF5 due to its limitation of natively storing complex number.

It's worth mention that [JLD](https://github.com/JuliaIO/JLD.jl) and
[JLD2](https://github.com/JuliaIO/JLD2.jl) are also excellent alternatives but they
storing the whole type information that may broke read when this package changes or
remove certain type definitions. Users may choose them as auxiliary tools.

## Public Interface
```@autodocs
Modules = [Quaycle]
Pages = ["h5getstore.jl", "h5solution.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```
