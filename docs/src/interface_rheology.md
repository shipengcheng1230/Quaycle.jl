# Rheology

This package implements plastic deformation as the key for modeling asthenosphere dynamics. Currently, only `DiffusionCreep` and `DislocationCreep` are supported.

## Public Interface
```@autodocs
Modules = [JuEQ]
Pages = ["rheology.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```

## References
Hirth, G., & Kohlstedt, D. (2003). Rheology of the Upper Mantle and the Mantle Wedge: A View from the Experimentalists. In Inside the Subduction Factory (pp. 83–105). American Geophysical Union (AGU). https://doi.org/10.1029/138GM06

Karato, S. (2010). Rheology of the Earth’s mantle: A historical review. Gondwana Research, 18(1), 17–45. https://doi.org/10.1016/j.gr.2010.03.004

Kohlstedt, D. L., & Hansen, L. N. (2015). 2.18 - Constitutive Equations, Rheological Behavior, and Viscosity of Rocks. In G. Schubert (Ed.), Treatise on Geophysics (Second Edition) (pp. 441–472). Oxford: Elsevier. https://doi.org/10.1016/B978-0-444-53802-4.00042-7
