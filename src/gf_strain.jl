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
function shear_traction_sbarbot_on_okada(::STRIKING, σvec::AbstractVector, dip::T) where T<:Real
    -σvec[2] * sind(dip) - σvec[5] * cosd(dip)
end

function shear_traction_sbarbot_on_okada(::DIPPING, σvec::AbstractVector, dip::T) where T<:Real
    (σvec[6] - σvec[1])/2 * sind(2dip) - σvec[3] * cosd(2dip)
end

function coordinate_sbarbot2okada!(u::AbstractVector)
    u[1], u[4] = u[4], u[1] # σxx, σyy = σ22, σ11
    u[3], u[5] = -u[5], -u[3] # σxz, σyz = -σ23, -σ13
end

"""
    stress_greens_func(ma::SBarbotMeshEntity{3}, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault,
        comp::NTuple{N, Symbol}; kwargs...) where {T, N}

Compute traction Green's function from [`SBarbotTet4MeshEntity`](@ref) or [`SBarbotHex8MeshEntity`](@ref) to [`RectOkadaMesh`](@ref)

## Arguments
- `ma::SBarbotMeshEntity{3}`: asthenosphere mesh
- `mf::RectOkadaMesh`: fault mesh
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ft::FlatPlaneFault`: fault type, either [`DIPPING()`](@ref) or [`STRIKING()`](@ref)
- `comp`: the strain component(s) to be considered. Either a singleton of `:xx`, `:xy`, `:xz`,
    `:yy`, `:yz`, `:zz` or tuple of a few ones

## Output
A tuple of ``n`` matrix, each represents interaction from one strain to the traction on fault.
"""
function stress_greens_func(ma::SBarbotMeshEntity{3}, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault, comp::NTuple{N, Symbol}; kwargs...) where {T, N}
    f = (c) -> stress_greens_func(ma, mf, λ, μ, ft, c; kwargs...)
    map(f, comp)
end

function stress_greens_func(ma::SBarbotMeshEntity{3}, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault, comp::Symbol; kwargs...) where T
    st = SharedArray{T}(mf.nx * mf.nξ, length(ma.tag))
    stress_greens_func!(st, ma, mf, λ, μ, ft, comp; kwargs...)
    return sdata(st)
end

function stress_greens_func_chunk!(
    st::SharedArray{T, 2}, subs::AbstractArray, ma::SBarbotMeshEntity{3}, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault, comp::Symbol;
    quadrature::Union{Nothing, NTuple}=nothing) where T

    ν = λ / 2 / (λ + μ)
    σ = Vector{T}(undef, 6)
    uϵ = unit_strain(Val(comp), T)
    i2s = CartesianIndices((mf.nx, mf.nξ))

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of fault, index of volume
        q = i2s[i]
        if isa(ma, SBarbotHex8MeshEntity)
            sbarbot_stress_hex8!(σ, mf.y[q[2]], mf.x[q[1]], -mf.z[q[2]], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
        elseif isa(ma, SBarbotTet4MeshEntity)
            sbarbot_stress_tet4!(σ, quadrature, mf.y[q[2]], mf.x[q[1]], -mf.z[q[2]], ma.A[j], ma.B[j], ma.C[j], ma.D[j], uϵ..., μ, ν)
        else
            error("Unsupported mesh entity type: $(typeof(ma)).")
        end
        st[i,j] = shear_traction_sbarbot_on_okada(ft, σ, mf.dip)
    end
end

"""
    stress_greens_func(ma::SBarbotMeshEntity{3}, λ::T, μ::T, comp::NTuple{N, Symbol};
        kwargs...) where {T, N}

Compute stress Green's function within [`SBarbotTet4MeshEntity`](@ref) or [`SBarbotHex8MeshEntity`](@ref)

## Arguments
- `ma::SBarbotMeshEntity{3}`: asthenosphere mesh
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `comp`: the strain component(s) to be considered. Either a singleton of `:xx`, `:xy`, `:xz`,
    `:yy`, `:yz`, `:zz` or tuple of a few ones

## Output
A tuple of ``n`` tuple of matrix, each tuple represents interaction from one strain, w.r.t. `comp`
    to the stress components, whose order is
    ``σ_{xx}``, ``σ_{xy}``, ``σ_{xz}``, ``σ_{yy}``, ``σ_{yz}``, ``σ_{zz}``,
    each of which is a matrix, within asthenosphere.
"""
function stress_greens_func(ma::SBarbotMeshEntity{3}, λ::T, μ::T, comp::NTuple{N, Symbol}; kwargs...) where {T, N}
    f = (c) -> stress_greens_func(ma, λ, μ, c; kwargs...)
    map(f, comp)
end

function stress_greens_func(ma::SBarbotMeshEntity{3}, λ::T, μ::T, comp::Symbol; kwargs...) where T
    numelements = length(ma.tag)
    st = ntuple(_ -> SharedArray{T}(numelements, numelements), Val(6))
    stress_greens_func!(st, ma, λ, μ, comp; kwargs...)
    return ntuple(x -> st[x] |> sdata, 6)
end

function stress_greens_func_chunk!(
    st::NTuple{N, <:SharedArray}, subs::AbstractArray, ma::SBarbotMeshEntity{3}, λ::T, μ::T, comp::Symbol;
    quadrature::Union{Nothing, NTuple}=nothing) where {T, N}

    ν = λ / 2 / (λ + μ)
    σ = Vector{T}(undef, 6)
    uϵ = unit_strain(Val(comp), T)
    indexST = Base.OneTo(6)

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of recv, index of src
        if isa(ma, SBarbotHex8MeshEntity)
            sbarbot_stress_hex8!(σ, ma.x1[i], ma.x2[i], ma.x3[i], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
        elseif isa(ma, SBarbotTet4MeshEntity)
            sbarbot_stress_tet4!(σ, quadrature, ma.x1[i], ma.x2[i], ma.x3[i], ma.A[j], ma.B[j], ma.C[j], ma.D[j], uϵ..., μ, ν)
        else
            error("Unsupported mesh entity type: $(typeof(ma)).")
        end
        coordinate_sbarbot2okada!(σ)
        for ic in indexST
            st[ic][i,j] = σ[ic]
        end
    end
end
