# Green's Function

The interactions among fault patches, asthenosphere elements as well as
  between fault and asthenosphere are computed via convolution of Green's
  Function. Currently supported Green's Function:
- 1D elastic line dislocation
- 2D elastic rectangular dislocation
- 2D elastic triangular dislocation
- 3D inelastic strain in Hex8 or Tet4 elements

Other types, such as 2D inelastic (plane stress or antiplane stress), polygon
dislocation, are WIP (**PR are welcome!**).

All functions translated here are ensured to be type stable
and have minimum allocation. Broadcasting isn't supported here
([why?](https://julialang.org/blog/2017/01/moredots)) so you can easily embed them
into multiprocessors parallel computation, which is implemented here.

Also be aware of the coordinate system difference among all the Green's function
  provided here. Users are encouraged to view the original sources for further
  details. This package preserves the original coordinate systems as well as
  function arguments for all, whose outcomes are also tested against original ones.

## Public Interface
```@autodocs
Modules = [Quaycle]
Pages = ["gf.jl", "gf_dislocation.jl", "gf_strain.jl", "gf_operator.jl",
         "okada_dc3d.jl", "sbarbot_hex8.jl", "sbarbot_tet4.jl", "nikkhoo_td.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```

## References
Okada, Y. (1992). Internal deformation due to shear and tensile faults in a half-space. Bulletin of the Seismological Society of America, 82(2), 1018–1040.

Pan, E., Yuan, J. H., Chen, W. Q., & Griffith, W. A. (2014). Elastic Deformation due to Polygonal Dislocations in a Transversely Isotropic Half‐SpaceElastic Deformation due to Polygonal Dislocations in a Transversely Isotropic Half‐Space. Bulletin of the Seismological Society of America, 104(6), 2698–2716. https://doi.org/10.1785/0120140161

Nikkhoo, M., & Walter, T. R. (2015). Triangular dislocation: an analytical, artefact-free solution. Geophysical Journal International, 201(2), 1119–1141. https://doi.org/10.1093/gji/ggv035

Barbot, S., Moore, J. D. P., & Lambert, V. (2017). Displacement and Stress Associated with Distributed Anelastic Deformation in a Half‐Space. Bulletin of the Seismological Society of America, 107(2), 821–855. https://doi.org/10.1785/0120160237

Barbot, S. (2018). Deformation of a Half‐Space from Anelastic Strain Confined in a Tetrahedral Volume. Bulletin of the Seismological Society of America, 108(5A), 2687–2712. https://doi.org/10.1785/0120180058
