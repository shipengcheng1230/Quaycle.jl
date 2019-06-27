# Green's Function

The interactions among fault patches, asthenosphere elements as well as
  between fault and asthenosphere are computed via convolution of Green's
  Function. Currently supported Green's Function:
- 1D/2D elastic dislocation
- 3D inelastic element of Hex8 or Tet4

Other types, such as 2D inelastic (plane stress or antiplane stress), curved
dislocation, are WIP.

## Public Interface
```@autodocs
Modules = [JuEQ]
Pages = ["gf.jl", "gf_okada.jl", "gf_sbarbot.jl", "gf_operator.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```

## References
Okada, Y. (1992). Internal deformation due to shear and tensile faults in a half-space. Bulletin of the Seismological Society of America, 82(2), 1018–1040.

Barbot, S., Moore, J. D. P., & Lambert, V. (2017). Displacement and Stress Associated with Distributed Anelastic Deformation in a Half‐Space. Bulletin of the Seismological Society of America, 107(2), 821–855. https://doi.org/10.1785/0120160237

Barbot, S. (2018). Deformation of a Half‐Space from Anelastic Strain Confined in a Tetrahedral Volume. Bulletin of the Seismological Society of America, 108(5A), 2687–2712. https://doi.org/10.1785/0120180058
