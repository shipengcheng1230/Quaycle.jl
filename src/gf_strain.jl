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

"Shear traction from output of SBarbot green's function on [`dc3d`](@ref)."
@inline function shear_traction_sbarbot_on_okada(::STRIKING, σvec::AbstractVector, dip::T) where T<:Real
    -σvec[2] * sind(dip) - σvec[5] * cosd(dip)
end

@inline function shear_traction_sbarbot_on_okada(::DIPPING, σvec::AbstractVector, dip::T) where T<:Real
    (σvec[6] - σvec[1])/2 * sind(2dip) - σvec[3] * cosd(2dip)
end

@inline function flip_stress_components_sbarbot!(u::AbstractVector)
    u[1], u[4] = u[4], u[1] # σxx, σyy = σ22, σ11
    u[3], u[5] = -u[5], -u[3] # σxz, σyz = -σ23, -σ13
end

"Shear traction from output of SBarbot on triangular mesh."
@inline function shear_traction_sbarbot_on_td(ft::FlatPlaneFault, σvec::V, ss::V, ds::V, ts::V) where V
    fx = σvec[4] * ts[1] + σvec[2] * ts[2] - σvec[5] * ts[3]
    fy = σvec[2] * ts[1] + σvec[1] * ts[2] - σvec[3] * ts[3]
    fz = -σvec[5] * ts[1] - σvec[3] * ts[2] + σvec[6] * ts[3]
    _shear_traction_td(ft, fx, fy, fz, ss, ds, ts)
end

"""
    stress_greens_func(ma::SBarbotMeshEntity{3}, mf::AbstractMesh{2}, λ::T, μ::T, ft::PlaneFault,
        ϵcomp::NTuple{N, Symbol}; kwargs...) where {T, N}

Compute traction Green's function from [`SBarbotTet4MeshEntity`](@ref) or [`SBarbotHex8MeshEntity`](@ref) to [`RectOkadaMesh`](@ref)

## Arguments
- `ma::SBarbotMeshEntity{3}`: asthenosphere mesh
- `mf::AbstractMesh{2}`: fault mesh
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ft::FlatPlaneFault`: fault type, either [`DIPPING()`](@ref) or [`STRIKING()`](@ref)
- `ϵcomp`: the strain ϵcomponent(s) to be considered

## Output
A tuple of ``n`` matrix, each represents interaction from one strain to the traction on fault.
"""
function stress_greens_func(ma::SBarbotMeshEntity{3}, mf::AbstractMesh{2}, λ::T, μ::T, ft::PlaneFault, ϵcomp::NTuple{N, Symbol}; kwargs...) where {T, N}
    f = (c) -> stress_greens_func(ma, mf, λ, μ, ft, c; kwargs...)
    map(f, ϵcomp)
end

function stress_greens_func(ma::SBarbotMeshEntity{3}, mf::AbstractMesh{2}, λ::T, μ::T, ft::PlaneFault, ϵcomp::Symbol; kwargs...) where T
    if isa(mf, RectOkadaMesh)
        num_patch = mf.nx * mf.nξ
    elseif isa(mf, TDTri3MeshEntity)
        num_patch = length(mf.tag)
    else
        error("Unsupported mesh entity type: $(typeof(mf)).")
    end
    st = SharedArray{T}(num_patch, length(ma.tag))
    stress_greens_func!(st, ma, mf, λ, μ, ft, ϵcomp; kwargs...)
    return sdata(st)
end

function stress_greens_func_chunk!(
    st::SharedArray{T, 2}, subs::AbstractArray, ma::SBarbotMeshEntity{3}, mf::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault, ϵcomp::Symbol;
    quadrature::Union{Nothing, NTuple}=nothing) where T

    ν = λ / 2 / (λ + μ)
    σ = Vector{T}(undef, 6)
    uϵ = unit_strain(Val(ϵcomp), T)
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

function stress_greens_func_chunk!(
    st::SharedArray{T, 2}, subs::AbstractArray, ma::SBarbotMeshEntity{3}, mf::TDTri3MeshEntity, λ::T, μ::T, ft::PlaneFault, ϵcomp::Symbol;
    quadrature::Union{Nothing, NTuple}=nothing) where T

    ν = λ / 2 / (λ + μ)
    σ = Vector{T}(undef, 6)
    uϵ = unit_strain(Val(ϵcomp), T)

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of fault, index of volume
        if isa(ma, SBarbotHex8MeshEntity)
            sbarbot_stress_hex8!(σ, mf.y[i], mf.x[i], -mf.z[i], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
        elseif isa(ma, SBarbotTet4MeshEntity)
            sbarbot_stress_tet4!(σ, quadrature, mf.y[i], mf.x[i], -mf.z[i], ma.A[j], ma.B[j], ma.C[j], ma.D[j], uϵ..., μ, ν)
        else
            error("Unsupported mesh entity type: $(typeof(ma)).")
        end
        st[i,j] = shear_traction_sbarbot_on_td(ft, σ, mf.ss[i], mf.ds[i], mf.ts[i])
    end
end

"""
    stress_greens_func(ma::SBarbotMeshEntity{3}, λ::T, μ::T,
        ϵcomp::NTuple{N1, Symbol}, σcomp::NTuple{N2, Symbol}; kwargs...) where {T, N1, N2}

Compute stress Green's function within [`SBarbotTet4MeshEntity`](@ref) or [`SBarbotHex8MeshEntity`](@ref)

## Arguments
- `ma::SBarbotMeshEntity{3}`: asthenosphere mesh
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ϵcomp`: the strain ϵcomponent(s) to be considered
- `σcomp::NTuple{N, Symbol}`: stress components to consider

## Output
A tuple of ``n`` tuple of matrix, each tuple represents interaction from one strain component, w.r.t. `ϵcomp`
    to the stress ϵcomponents, whose order is ``σ_{ij}`` whose order is the same as `σcomp`,
    each of which is a matrix.
"""
function stress_greens_func(ma::SBarbotMeshEntity{3}, λ::T, μ::T, ϵcomp::NTuple{N1, Symbol}, σcomp::NTuple{N2, Symbol}; kwargs...) where {T, N1, N2}
    f = (x) -> stress_greens_func(ma, λ, μ, x, σcomp; kwargs...)
    map(f, ϵcomp)
end

function stress_greens_func(ma::SBarbotMeshEntity{3}, λ::T, μ::T, ϵcomp::Symbol, σcomp::NTuple{N, Symbol}; kwargs...) where {T, N}
    numelements = length(ma.tag)
    numσ = length(σcomp)
    st = ntuple(_ -> SharedArray{T}(numelements, numelements), Val(numσ))
    stress_greens_func!(st, ma, λ, μ, ϵcomp, σcomp; kwargs...)
    return ntuple(x -> st[x] |> sdata, numσ)
end

function stress_greens_func_chunk!(
    st::NTuple{N, <:SharedArray}, subs::AbstractArray, ma::SBarbotMeshEntity{3}, λ::T, μ::T, ϵcomp::Symbol, σcomp::NTuple{N, Symbol};
    quadrature::Union{Nothing, NTuple}=nothing) where {T, N}

    ν = λ / 2 / (λ + μ)
    σ = Vector{T}(undef, 6)
    uϵ = unit_strain(Val(ϵcomp), T)
    σindex = [_ϵσindex(Val(x)) for x in σcomp]
    indexST = Base.OneTo(length(σcomp))

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of recv, index of src
        if isa(ma, SBarbotHex8MeshEntity)
            sbarbot_stress_hex8!(σ, ma.x1[i], ma.x2[i], ma.x3[i], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
        elseif isa(ma, SBarbotTet4MeshEntity)
            sbarbot_stress_tet4!(σ, quadrature, ma.x1[i], ma.x2[i], ma.x3[i], ma.A[j], ma.B[j], ma.C[j], ma.D[j], uϵ..., μ, ν)
        else
            error("Unsupported mesh entity type: $(typeof(ma)).")
        end
        flip_stress_components_sbarbot!(σ)
        for ic in indexST
            st[ic][i,j] = σ[σindex[ic]]
        end
    end
end

## displacement - strain Green's function
function disp_greens_func(x::T, y::T, z::T, mf::SBarbotMeshEntity, ν::T) where T<:Real
    # TODO
end
