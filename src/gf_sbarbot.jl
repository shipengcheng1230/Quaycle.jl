"""
Given the axes transformation rule that ``x ⟶ x_2, \\; y ⟶ x_1, \\; z ⟶ -x_3``,
the corresponding strain mapping is:
```math
ϵ_{xx} ⟶ ϵ_{22}
```
```math
ϵ_{xy} ⟶ ϵ_{12}
```
```math
ϵ_{xz} ⟶ -ϵ_{23}
```
```math
ϵ_{yy} ⟶ ϵ_{11}
```
```math
ϵ_{yz} ⟶ -ϵ_{13}
```
```math
ϵ_{zz} ⟶ ϵ_{33}
```
"""
@inline unit_strain(::Val{:xx}, T=Float64) = [zero(T), zero(T), zero(T), one(T), zero(T), zero(T)]
@inline unit_strain(::Val{:xy}, T=Float64) = [zero(T), one(T), zero(T), zero(T), zero(T), zero(T)]
@inline unit_strain(::Val{:xz}, T=Float64) = [zero(T), zero(T), zero(T), zero(T), -one(T), zero(T)]
@inline unit_strain(::Val{:yy}, T=Float64) = [one(T), zero(T), zero(T), zero(T), zero(T), zero(T)]
@inline unit_strain(::Val{:yz}, T=Float64) = [zero(T), zero(T), -one(T), zero(T), zero(T), zero(T)]
@inline unit_strain(::Val{:zz}, T=Float64) = [zero(T), zero(T), zero(T), zero(T), zero(T), one(T)]

"Shear traction from output of SBarbot green's function."
function shear_traction_sbarbot(::STRIKING, σvec::AbstractVector, λ::T, μ::T, dip::T) where T<:Real
    σvec[2] * sind(dip) + σvec[5] * cosd(dip)
end

function shear_traction_sbarbot(::DIPPING, σvec::AbstractVector, λ::T, μ::T, dip::T) where T<:Real
    (σvec[1] - σvec[6])/2 * sind(2dip) + σvec[3] * cosd(2dip)
end
