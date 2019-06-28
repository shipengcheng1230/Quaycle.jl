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
- `mesh::RectOkadaMesh`: the line mesh coupled with [`dc3d`](@ref)
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

    function __convert_to_fourier_domain__()
        FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
        x1 = zeros(T, 2 * mesh.nx - 1)
        p1 = plan_rfft(x1, flags=parameters["FFT"]["FLAG"])
        st_dft = Array{Complex{T}}(undef, mesh.nx, mesh.nξ, mesh.nξ)
        @inbounds for l = 1: mesh.nξ, j = 1: mesh.nξ
            # reference: https://github.com/JuliaMatrices/ToeplitzMatrices.jl/blob/cbe29c344be8363f33eb17090121f8cff600b72e/src/ToeplitzMatrices.jl#L627
            st_dft[:,j,l] .= p1 * [st[:,j,l]; reverse(st[2:end,j,l])]
        end
        return st_dft
    end

    fourier_domain ? __convert_to_fourier_domain__() : sdata(st)
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

function stress_greens_func_chunk!(st::SharedArray{T, 3}, subs::AbstractArray, mesh::RectOkadaMesh, λ::T, μ::T, ft::FlatPlaneFault; nrept::Integer=2, buffer_ratio::Real=zero(T)) where T
    ud = unit_dislocation(ft, T)
    lrept = (buffer_ratio + one(T)) * (mesh.Δx * mesh.nx)
    u = Vector{T}(undef, 12)
    α = (λ + μ) / (λ + 2μ)
    @inbounds for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        okada_gf_periodic_bc!(u, mesh.x[i], mesh.y[j], mesh.z[j], α, mesh.dep, mesh.dip, mesh.ax[1], mesh.aξ[l], ud, nrept, lrept)
        st[i,j,l] = shear_traction_dc3d(ft, u, λ, μ, mesh.dip)
    end
end

"""
Periodic summation of green's function.
- `lrept` represents jumping periodic block for simulating locked region at along-strike edge.
- `nrept` represents number of blocks to be summed.
"""
function okada_gf_periodic_bc!(u::AbstractVector{T}, x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    nrept::Integer, lrept::T) where {T<:Real}
    fill!(u, zero(T))
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d(x, y, z, α, depth, dip, al .+ i * lrept, aw, disl)
    end
end

"Dislocation direction: *hanging wall* - *foot wall*."
@inline unit_dislocation(::DIPPING, T=Float64) = [zero(T), one(T), zero(T)]
@inline unit_dislocation(::STRIKING, T=Float64) = [one(T), zero(T), zero(T)]

"Normal of hanging outwards: ``(0,\\; -\\sin{θ},\\; \\cos{θ})``."
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

## okada displacement - stress green's function, src on the fault but recv in the other volume
"Compute 6 stress components from [`dc3d`](@ref) output."
@inline function stress_components!(σ::T, u::T, λ::U, μ::U) where {T<:AbstractVector, U<:Real}
    ϵkk = u[4] + u[8] + u[12]
    σ[1] = λ * ϵkk + 2μ * u[4]
    σ[2] = μ * (u[5] + u[7])
    σ[3] = μ * (u[6] + u[10])
    σ[4] = λ * ϵkk + 2μ * u[8]
    σ[5] = μ * (u[9] + u[11])
    σ[6] = λ * ϵkk + 2μ * u[12]
end

function stress_components(u::T, λ::U, μ::U) where {T<:AbstractVector, U<:Real}
    σ = Vector{eltype(u)}(undef, 6)
    stress_components!(σ, u, λ, μ)
    return σ
end

"""
    stress_greens_func(mf::RectOkadaMesh, ma::SBarbotMeshEntity{3}, λ::T, μ::T, ft::FlatPlaneFault;
        kwargs...) where {T<:Real, I<:Integer}

Compute stress Green's function from [`RectOkadaMesh`](@ref) to [`SBarbotTet4MeshEntity`](@ref) or [`SBarbotHex8MeshEntity`](@ref)

## Arguments
- `mf::RectOkadaMesh`: mesh of fault
- `ma::SBarbotMeshEntity{3}`: mesh of asthenosphere
- `λ::T`: Lamé's first parameter
- `μ::T`: shear modulus
- `ft::FlatPlaneFault`: fault type, either [`DIPPING()`](@ref) or [`STRIKING()`](@ref)

## Output
The output is a tuple of 6 matrix, each corresponds ``σ_{xx}``, ``σ_{xy}``, ``σ_{xz}``,
    ``σ_{yy}``, ``σ_{yz}``, ``σ_{zz}``
"""
function stress_greens_func(mf::RectOkadaMesh, ma::SBarbotMeshEntity{3}, λ::T, μ::T, ft::FlatPlaneFault; kwargs...) where {T<:Real, I<:Integer}
    st = ntuple(_ -> SharedArray{T}(length(ma.tag), mf.nx * mf.nξ), Val(6))
    stress_greens_func!(st, mf, ma, λ, μ, ft; kwargs...)
    return ntuple(x -> st[x] |> sdata, 6)
end

function stress_greens_func_chunk!(
    st::NTuple{N, <:SharedArray}, subs::AbstractArray, mf::RectOkadaMesh, ma::SBarbotMeshEntity{3},
    λ::T, μ::T, ft::FlatPlaneFault; nrept::Integer=2, buffer_ratio::Real=0.0,
    ) where {T<:Real, N, I<:Integer}
    ud = unit_dislocation(ft)
    lrept = (buffer_ratio + one(T)) * (mf.Δx * mf.nx)
    α = (λ + μ) / (λ + 2μ)
    u = Vector{T}(undef, 12)
    σ = Vector{T}(undef, 6)
    i2s = CartesianIndices((mf.nx, mf.nξ))
    indexST = Base.OneTo(6)

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of volume, index of fault
        q = i2s[j] # return (ix, iξ)
        okada_gf_periodic_bc!(u, ma.x2[i], ma.x1[i], -ma.x3[i], α, mf.dep, mf.dip, mf.ax[q[1]], mf.aξ[q[2]], ud, nrept, lrept)
        stress_components!(σ, u, λ, μ)
        for ind in indexST
            st[ind][i,j] = σ[ind]
        end
    end
end
