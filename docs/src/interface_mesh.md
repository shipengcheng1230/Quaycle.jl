# Mesh

This package provides some buildin structured mesh functionality and some
  data structures coupled with existing Green's functions where users could
  explore external mesh tools.

This package also has rich utilities based on [Gmsh](http://gmsh.info/).
  To use them, it is recommended to install
  [GmshTools.jl](https://github.com/shipengcheng1230/GmshTools.jl) which
  automatically downloads Gmsh SDK and provides a few convenient macros.


```julia
pkg> add GmshTools
julia> using JuEQ
julia> using GmshTools
```

Users are encouraged to read [Gmsh Julia API](https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.jl) for more comprehensive operations.

## Public Interface
```@autodocs
Modules = [JuEQ]
Pages = ["mesh.jl", "gmshtools.jl", "space.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```

## References
Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: A 3-D finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering, 79(11), 1309–1331. https://doi.org/10.1002/nme.2579
