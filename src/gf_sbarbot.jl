export sbarbot_stress_gf_tensor

@gen_shared_chunk_call sbarbot_stress_gf_tensor

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

function sbarbot_stress_gf_tensor(ma::SBarbotHex8MeshEntity, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault, comp::NTuple{N, Symbol}) where {T, N}
    f = (c) -> sbarbot_stress_gf_tensor(ma, mf, λ, μ, ft, c)
    map(f, comp)
end

function sbarbot_stress_gf_tensor(ma::SBarbotHex8MeshEntity, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault, comp::Symbol) where T
    st = SharedArray{T}(mf.nx * mf.nξ, length(ma.tag))
    sbarbot_stress_gf_tensor!(st, ma, mf, λ, μ, ft, comp)
    return st
end

function sbarbot_stress_gf_tensor_chunk!(
    st::SharedArray{T, 2}, subs::AbstractArray, ma::SBarbotHex8MeshEntity, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault, comp::Symbol,
    ) where T

    ν = λ / 2 / (λ + μ)
    σ = Vector{T}(undef, 6)
    uϵ = unit_strain(Val(comp), T)
    i2s = CartesianIndices((mf.nx, mf.nξ))

    for sub in subs
        i, j = sub[1], sub[2] # index of fault, index of volume
        q = i2s[i]
        sbarbot_stress_hex8!(σ, mf.y[q[2]], mf.x[q[1]], -mf.z[q[2]], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
        st[i,j] = shear_traction_sbarbot(ft, σ, λ, μ, mf.dip)
    end
end
