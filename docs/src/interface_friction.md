# Friction

This package adopts what is called **rate-and-state** friction as one of the
essential components for modeling fault dynamics. Currently, we only support
single state variable ``θ``.

## Public Interface
```@autodocs
Modules = [Quaycle]
Pages = ["friction.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```

## References
Dieterich, J. (1979). Modeling of rock friction: 1. Experimental results and constitutive equations. Journal of Geophysical Research: Solid Earth, 84(B5), 2161–2168. https://doi.org/10.1029/JB084iB05p02161

Ruina, A. (1983). Slip instability and state variable friction laws. Journal of Geophysical Research: Solid Earth, 88(B12), 10359–10370. https://doi.org/10.1029/JB088iB12p10359

Rubin, A. M., & Ampuero, J.-P. (2005). Earthquake nucleation on (aging) rate and state faults. Journal of Geophysical Research: Solid Earth, 110(B11). https://doi.org/10.1029/2005JB003686
