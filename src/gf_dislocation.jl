## okada displacement - traction green's function, src & recv both on the same fault mesh
"""
    stress_greens_func(mesh::LineOkadaMesh, λ::T, μ::T, ft::FlatPlaneFault; kwargs...) where T

Compute traction Green function in 1-D elastic fault in [`LineOkadaMesh`](@ref).

## Arguments
- `mesh::LineOkadaMesh`: the line mesh coupled with [`dc3d`](@ref)
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ft::FlatPlaneFault`: fault type, either [`DIPPING()`](@ref) or [`STRIKING()`](@ref)

### KWARGS Arguments
- `ax_ratio::Real`: ratio of along-strike to along-downdip, default is `12.5`
"""
function stress_greens_func(mesh::LineOkadaMesh, λ::T, μ::T, ft::FlatPlaneFault; kwargs...) where T
    st = SharedArray{T}(mesh.nξ, mesh.nξ)
    stress_greens_func!(st, mesh, λ, μ, ft; kwargs...)
    return sdata(st)
end

"""
    stress_greens_func(mesh::RectOkadaMesh, λ::T, μ::T, ft::FlatPlaneFault;
        fourier_domain=true, kwargs...) where T

Compute traction Green's function in 2-D elastic fault in [`RectOkadaMesh`](@ref). Translational symmetry is considered.

## Arguments
- `mesh::RectOkadaMesh`: the rectangular mesh coupled with [`dc3d`](@ref)
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ft::FlatPlaneFault`: fault type, either [`DIPPING()`](@ref) or [`STRIKING()`](@ref)

### KWARGS Arguments
- `fourier_domain::Bool`: whether or not transform the tensor to fourier domain
- `nrept::Integer`: number of periodic summation is performed
- `buffer_ratio::Real`: ratio of length of buffer zone (along-strike) to that of fault (along-strike)
    It is recommended to set at least `1` for strike-slip fault for mimicing zero-dislocation at ridge on both sides.
"""
function stress_greens_func(mesh::RectOkadaMesh, λ::T, μ::T, ft::FlatPlaneFault; fourier_domain=true, kwargs...) where T
    st = SharedArray{T}(mesh.nx, mesh.nξ, mesh.nξ)
    stress_greens_func!(st, mesh, λ, μ, ft; kwargs...)

    if fourier_domain
        x1 = zeros(T, 2 * mesh.nx - 1)
        p1 = plan_rfft(x1, flags=parameters["FFT"]["FLAG"])
        st_dft = Array{Complex{T}}(undef, mesh.nx, mesh.nξ, mesh.nξ)
        @inbounds for l = 1: mesh.nξ, j = 1: mesh.nξ
            # reference: https://github.com/JuliaMatrices/ToeplitzMatrices.jl/blob/cbe29c344be8363f33eb17090121f8cff600b72e/src/ToeplitzMatrices.jl#L627
            st_dft[:,j,l] .= p1 * [st[:,j,l]; reverse(st[2:end,j,l])]
        end
        return st_dft
    else
        sdata(st)
    end
end

function stress_greens_func_chunk!(st::SharedArray{T, 2}, subs::AbstractArray, mesh::LineOkadaMesh, λ::T, μ::T, ft::FlatPlaneFault; ax_ratio::Real=12.5) where T
    ud = unit_dislocation(ft, T)
    ax = mesh.nξ * mesh.Δξ * ax_ratio .* [-one(T), one(T)]
    α = (λ + μ) / (λ + 2μ)
    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2]
        u = dc3d(mesh.x, mesh.y[i], mesh.z[i], α, mesh.dep, mesh.dip, ax, mesh.aξ[j], ud)
        st[i,j] = shear_traction_dc3d(ft, u, λ, μ, mesh.dip)
    end
end

function stress_greens_func_chunk!(st::SharedArray{T, 3}, subs::AbstractArray, mesh::RectOkadaMesh, λ::T, μ::T, ft::FlatPlaneFault; nrept::Integer=2, buffer_ratio::Real=0) where T
    @assert buffer_ratio ≥ 0 "Argument `buffer_ratio` must be ≥ 0."
    ud = unit_dislocation(ft, T)
    lrept = (buffer_ratio + one(T)) * (mesh.Δx * mesh.nx)
    u = Vector{T}(undef, 12)
    α = (λ + μ) / (λ + 2μ)
    @inbounds for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        dc3d_periodic_bc!(u, mesh.x[i], mesh.y[j], mesh.z[j], α, mesh.dep, mesh.dip, mesh.ax[1], mesh.aξ[l], ud, nrept, lrept)
        st[i,j,l] = shear_traction_dc3d(ft, u, λ, μ, mesh.dip)
    end
end

"""
Periodic summation of green's function.
- `lrept` represents jumping periodic block for simulating locked region at along-strike edge.
- `nrept` represents number of blocks to be summed.
"""
function dc3d_periodic_bc!(u::AbstractVector{T}, x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    nrept::Integer, lrept::T) where {T<:Real}
    fill!(u, zero(T))
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d(x, y, z, α, depth, dip, al .+ i * lrept, aw, disl)
    end
end

"Dislocation direction: *hanging wall* - *foot wall*."
@inline unit_dislocation(::DIPPING, T=Float64) = [zero(T), one(T), zero(T)]
@inline unit_dislocation(::STRIKING, T=Float64) = [one(T), zero(T), zero(T)]

"Normal of hanging inwards: ``(0,\\; -\\sin{θ},\\; \\cos{θ})``."
@inline function shear_traction_dc3d(::DIPPING, u::AbstractVector, λ::T, μ::T, dip::T) where T
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    (σzz - σyy)/2 * sind(2dip) + τyz * cosd(2dip)
end

@inline function shear_traction_dc3d(::STRIKING, u::AbstractVector, λ::T, μ::T, dip::T) where T
    σxy = μ * (u[5] + u[7])
    σxz = μ * (u[6] + u[10])
    -σxy * sind(dip) + σxz * cosd(dip)
end

"""
    stress_greens_func(mesh::TDTri3MeshEntity, λ::T, μ::T, ft::FlatPlaneFault; kwargs...) where T

Compute traction Green's function in 2-D elastic fault in [`TDTri3MeshEntity`](@ref).

## Arguments
- `mesh::TDTri3MeshEntity`: the triangular mesh coupled with [`td_stress_hs`](@ref)
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ft::FlatPlaneFault`: fault type, either [`DIPPING()`](@ref) or [`STRIKING()`](@ref)

### KWARGS Arguments
- `nrept::Integer`: number of periodic summation is performed
- `buffer_ratio::Real`: ratio of length of buffer zone (along-strike) to that of fault (along-strike).
    Notice the direction of strike is the average value of `mesh.ss` and the length is the strike projected
    maximum horizontal expansion.
"""
function stress_greens_func(mesh::TDTri3MeshEntity, λ::T, μ::T, ft::FlatPlaneFault; kwargs...) where T
    nume = length(mesh.tag)
    st = SharedArray{T}(nume, nume)
    stress_greens_func!(st, mesh, λ, μ, ft; kwargs...)
    return sdata(st)
end

function stress_greens_func_chunk!(st::SharedArray{T, 2}, subs::AbstractArray, mesh::TDTri3MeshEntity,
    λ::T, μ::T, ft::FlatPlaneFault; nrept::Integer=0, buffer_ratio::Real=0) where T

    vrept = nrept == 0 ? zeros(T, 3) : td_periodic_vec(mesh, buffer_ratio)
    ud = unit_dislocation(ft, T)
    σvec = Vector{T}(undef, 6)
    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of recv, index of src
        td_periodic_bc!(σvec, mesh.x[i], mesh.y[i], mesh.z[i], mesh.A[j], mesh.B[j], mesh.C[j], ud..., λ, μ, nrept, vrept)
        st[i,j] = shear_traction_td(ft, σvec, mesh.ss[i], mesh.ds[i], mesh.ts[i])
    end
end

@inline _shear_traction_td(::DIPPING, fx::T, fy::T, fz::T, ss::V, ds::V, ts::V) where {T, V} = fx * ds[1] + fy * ds[2] + fz * ds[3]
@inline _shear_traction_td(::STRIKING, fx::T, fy::T, fz::T, ss::V, ds::V, ts::V) where {T, V} = fx * ss[1] + fy * ss[2] + fz * ss[3]

"Compute shear traction from Tri3 dislocation."
@inline function shear_traction_td(ft::FlatPlaneFault, σvec::AbstractVector, ss::V, ds::V, ts::V) where {T, V}
    fx = σvec[1] * ts[1] + σvec[4] * ts[2] + σvec[5] * ts[3]
    fy = σvec[4] * ts[1] + σvec[2] * ts[2] + σvec[6] * ts[3]
    fz = σvec[5] * ts[1] + σvec[6] * ts[2] + σvec[3] * ts[3]
    _shear_traction_td(ft, fx, fy, fz, ss, ds, ts)
end

## okada displacement - stress green's function, src on the fault but recv in the other volume
"Compute 6 stress components from [`dc3d`](@ref) output."
@inline function stress_components_dc3d!(σ::T, u::T, λ::U, μ::U) where {T<:AbstractVector, U<:Real}
    ϵkk = u[4] + u[8] + u[12]
    σ[1] = λ * ϵkk + 2μ * u[4] # σxx
    σ[2] = μ * (u[5] + u[7]) # σxy
    σ[3] = μ * (u[6] + u[10]) # σxz
    σ[4] = λ * ϵkk + 2μ * u[8] # σyy
    σ[5] = μ * (u[9] + u[11]) # σyz
    σ[6] = λ * ϵkk + 2μ * u[12] # σzz
end

function stress_components_dc3d(u::T, λ::U, μ::U) where {T<:AbstractVector, U<:Real}
    σ = Vector{eltype(u)}(undef, 6)
    stress_components_dc3d!(σ, u, λ, μ)
    return σ
end

"Flip the output from [`td_stress_hs`](@ref) in accordance to [`stress_components_dc3d`](@ref)."
@inline function flip_stress_components_td!(σ::T) where T<:AbstractVector
    # σxx, σyy, σzz, σxy, σxz, σyz ⟶ σxx, σxy, σxz, σyy, σyz, σzz
    σ[2], σ[4] = σ[4], σ[2]
    σ[3], σ[6] = σ[6], σ[3]
    σ[3], σ[5] = σ[5], σ[3]
end

function td_periodic_bc!(σ::V,
    X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T,
    nrept::I, vrept::V) where {T, V, I}

    fill!(σ, zero(T))
    @fastmath @simd for i ∈ -nrept: nrept
        σ .+= td_stress_hs(X, Y, Z, P1 + i * vrept, P2 + i * vrept, P3 + i * vrept, Ss, Ds, Ts, λ, μ)
    end
end

function td_periodic_vec(mt::TDTri3MeshEntity, buffer_ratio::Real=0)
    @assert buffer_ratio ≥ 0 "Argument `buffer_ratio` must be ≥ 0."
    # maximum horizontal expansion project on the average strike direction
    rangeX = map(c -> extrema(map(x -> x[1], c)) |> y -> y[2] - y[1], [mt.A, mt.B, mt.C]) |> maximum
    rangeY = map(c -> extrema(map(x -> x[2], c)) |> y -> y[2] - y[1], [mt.A, mt.B, mt.C]) |> maximum
    avgss = reduce(+, mt.ss) |> normalize
    vrept = [rangeX * avgss[1], rangeY * avgss[2], zero(eltype(avgss))] * (buffer_ratio + 1)
end

"""
    stress_greens_func(mf::AbstractMesh{2}, ma::SBarbotMeshEntity{3},
        λ::T, μ::T, ft::FlatPlaneFault,
        σcomp::NTuple{N, Symbol}; kwargs...) where {T<:Real, I<:Integer, N}

Compute stress Green's function from fault mesh to asthenosphere mesh.

## Arguments
- `mf::AbstractMesh{2}`: mesh of fault
- `ma::SBarbotMeshEntity{3}`: mesh of asthenosphere
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ft::FlatPlaneFault`: fault type, either [`DIPPING()`](@ref) or [`STRIKING()`](@ref)
- `σcomp::NTuple{N, Symbol}`: stress components to consider

### KWARGS Arguments
The same as previously mentioned:
- `nrept::Integer`
- `buffer_ratio::Real`

## Output
The output is a tuple of `length(σcomp)` matrix, each corresponds ``σ_{ij}`` in the same order
    as given by `σcomp`.
"""
function stress_greens_func(mf::AbstractMesh{2}, ma::SBarbotMeshEntity{3}, λ::T, μ::T, ft::FlatPlaneFault, σcomp::NTuple{N, Symbol}; kwargs...) where {T<:Real, I<:Integer, N}
    if isa(mf, RectOkadaMesh)
        num_patch = mf.nx * mf.nξ
    elseif isa(mf, TDTri3MeshEntity)
        num_patch = length(mf.tag)
    else
        error("Unsupported mesh entity type: $(typeof(mf)).")
    end
    numσ = length(σcomp)
    st = ntuple(_ -> SharedArray{T}(length(ma.tag), num_patch), Val(numσ))
    stress_greens_func!(st, mf, ma, λ, μ, ft, σcomp; kwargs...)
    return map(sdata, st)
end

function stress_greens_func_chunk!(
    st::NTuple{N, <:SharedArray}, subs::AbstractArray, mf::RectOkadaMesh, ma::SBarbotMeshEntity{3},
    λ::T, μ::T, ft::FlatPlaneFault, σcomp::NTuple{N, Symbol}; nrept::Integer=2, buffer_ratio::Real=0.0) where {T<:Real, N, I<:Integer}

    ud = unit_dislocation(ft)
    lrept = (buffer_ratio + one(T)) * (mf.Δx * mf.nx)
    α = (λ + μ) / (λ + 2μ)
    u = Vector{T}(undef, 12)
    σ = Vector{T}(undef, 6)
    σindex = [_ϵσindex(Val(x)) for x in σcomp]
    i2s = CartesianIndices((mf.nx, mf.nξ))
    indexST = Base.OneTo(length(σcomp))

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of volume, index of fault
        q = i2s[j] # return (ix, iξ)
        dc3d_periodic_bc!(u, ma.x2[i], ma.x1[i], -ma.x3[i], α, mf.dep, mf.dip, mf.ax[q[1]], mf.aξ[q[2]], ud, nrept, lrept)
        stress_components_dc3d!(σ, u, λ, μ)
        for ind in indexST
            st[ind][i,j] = σ[σindex[ind]]
        end
    end
end

function stress_greens_func_chunk!(
    st::NTuple{N, <:SharedArray}, subs::AbstractArray, mf::TDTri3MeshEntity, ma::SBarbotMeshEntity{3},
    λ::T, μ::T, ft::FlatPlaneFault, σcomp::NTuple{N, Symbol}; nrept::Integer=0, buffer_ratio::Real=0) where {T<:Real, N, I<:Integer}

    vrept = nrept == 0 ? zeros(T, 3) : td_periodic_vec(mf, buffer_ratio)
    ud = unit_dislocation(ft)
    σ = Vector{T}(undef, 6)
    σindex = [_ϵσindex(Val(x)) for x in σcomp]
    indexST = Base.OneTo(length(σcomp))

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of volume, index of fault
        td_periodic_bc!(σ, ma.x2[i], ma.x1[i], -ma.x3[i], mf.A[j], mf.B[j], mf.C[j], ud..., λ, μ, nrept, vrept)
        flip_stress_components_td!(σ)
        for ind in indexST
            st[ind][i,j] = σ[σindex[ind]]
        end
    end
end

## displacement - displacement Green's function
function disp_greens_func(x::T, y::T, z::T, mf::RectOkadaMesh, α::T; nrept::Integer=2, buffer_ratio::Real=0.0) where T<:Real
    # TODO
end

function disp_greens_func(x::T, y::T, z::T, mf::LineOkadaMesh, α::T; ax_ratio::Real=12.5) where T<:Real
    # TODO
end

function disp_greens_func(x::T, y::T, z::T, mf::TDTri3MeshEntity, α::T; nrept::Integer=2, buffer_ratio::Real=0.0) where T<:Real
    # TODO
end
