## displacement traits definition
@traitdef IsDisplacementOnYZPlane{faulttype}
@traitimpl IsDisplacementOnYZPlane{faulttype} <- is_diplacement_on_yz_plane(faulttype)

is_diplacement_on_yz_plane(::Type{NormalFault}) = true
is_diplacement_on_yz_plane(::Type{ThrustFault}) = true
is_diplacement_on_yz_plane(::Type{StrikeSlipFault}) = false


## Grids and discretizing methods
abstract type AbstractDifferenceGrids{dim} end
abstract type BoundaryElementGrid{dim} <: AbstractDifferenceGrids{dim} end

struct BEMGrid_1D{U1<:AbstractVector, U2<:AbstractArray{<:AbstractVector}, T<:Number, I<:Integer} <: BoundaryElementGrid{1}
    ξ::U1 # along-downdip
    Δξ::T
    nξ::I
    aξ::U2 # cell boundary along-downdip (width)
    y::U1 # y-component of `ξ`
    z::U1 # z-component of `ξ`
    ax::U1 # (psudo) cell boundary along-strike (length)
end

struct BEMGrid_2D{U1<:AbstractVector, U2<:AbstractArray{<:AbstractVector}, T<:Number, I<:Integer} <: BoundaryElementGrid{2}
    x::U1 # along-strike
    ξ::U1 # along-downdip
    Δx::T
    Δξ::T
    nx::I
    nξ::I
    ax::U2 # cell boundary along-strike (length)
    aξ::U2 # cell boundary along-downdip (width)
    y::U1 # y-component of `ξ`
    z::U1 # z-component of `ξ`
    buffer_ratio::Union{T,I} # ratio of buffer zone length against along-strike length
end

"""
    discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ::T; ax_ratio=12.5)

Generate the grid for given 1D fault domain. The grids will be forced to start at (x=0, y=0, z=0).

## Arguments
- `Δξ`: grid space along-downdip
- `ax_ratio::Number`: ration of along-strike length agsinst along-downdip length for mimicing an extended
    2d (x & ξ) fault represented by 1d (ξ) domain. Default `ax_ratio=12.5` is more than enough for producing consistent results.
"""
function discretize(fa::PlaneFaultDomain{ftype, 1, T}, Δξ::T; ax_ratio=12.5) where {ftype<:PlaneFault, T<:Number}
    ξ, nξ, aξ = __divide_segment__(Val(:halfspace), promote(fa[:ξ], Δξ)...)
    ax = fa[:ξ] * ax_ratio .* [-one(T), one(T)]
    BEMGrid_1D(ξ, Δξ, nξ, aξ, ξ.*cospi(fa.dip/180), ξ.*sinpi(fa.dip/180), ax)
end

"""
    discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ; buffer=:auto) where {ftype <: PlaneFault}

Generate the grid for given 2D fault domain. The grids will be forced to start at (z=0) and spread symmetrically along x-axis w.r.t y-z plane.
    By such setting, we would be able to utilize the symmetry properties of stiffness tensor for performance speed up.

## Arguments
- `Δx, Δξ`: grid space along-strike and along-downdip respectively
- `buffer_ratio::Number: ration of buffer size against along-strike length for introducing *zero-dislocation* area at along-strike edges of defined fault domain.
"""
function discretize(fa::PlaneFaultDomain{ftype, 2, T}, Δx::T, Δξ::T; buffer_ratio=0.0::Number) where {ftype<:PlaneFault, T<:Number}

    ξ, nξ, aξ = __divide_segment__(Val(:halfspace), promote(fa[:ξ], Δξ)...)
    x, nx, ax = __divide_segment__(Val(:symmatzero), promote(fa[:x], Δx)...)

    BEMGrid_2D(x, ξ, Δx, Δξ, nx, nξ, ax, aξ, ξ.*cospi(fa.dip/180), ξ.*sinpi(fa.dip/180), buffer_ratio)
end

discretize(fa::PlaneFaultDomain{ftype, 1}; nξ::Integer=500, kwargs...) where {ftype<:PlaneFault} = discretize(fa, fa[:ξ]/nξ; kwargs...)
discretize(fa::PlaneFaultDomain{ftype, 2}; nx::Integer=60, nξ::Integer=30, kwargs...) where  {ftype<:PlaneFault} = discretize(fa, fa[:x]/nx, fa[:ξ]/nξ; kwargs...)

function __divide_segment__(::Val{:halfspace}, x::T, Δx::T) where {T<:Number}
    xi = range(zero(T), stop=-x+Δx, step=-Δx) .- Δx/2 |> collect
    xs = __segment_boarder__(xi, Δx)
    return xi, length(xi), xs
end

function __divide_segment__(::Val{:symmatzero}, x::T, Δx::T) where {T<:Number}
    xi = range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx) |> collect
    xs = __segment_boarder__(xi, Δx)
    return xi, length(xi), xs
end

__segment_boarder__(xi::AbstractVector{<:Number}, Δx::Number) = [[xx - Δx/2, xx + Δx/2] for xx in xi]


## Phyiscal properties structs
abstract type AbstractElasticProperties end

"""
Okada's [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) only applies on isotropic materials,
    therefore, elastic modulus are constrained to be scalars.
"""
@with_kw struct HomogeneousElasticProperties{T<:Number} <: AbstractElasticProperties
    λ::T # Lamé's first parameter
    μ::T # shear modulus
    α::T = (λ + μ) / (λ + 2μ) # material constants which equals (λ + μ) / (λ + 2μ)
end

@with_kw struct PlaneMaterialProperties{N, T<:Number, U<:AbstractVecOrMat, P<:AbstractArray} <: AbstractMaterialProperties{N}
    dims::Dims{N}
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    k::P # stiffness tensor
    σ::U # effective normal stress
    η::U # radiation damping
    vpl::T # plate rate, unlike pure Rate-State Friction simulation, here is restrained to be constant
    f0::T # ref. frictional coeff
    v0::T # ref. velocity

    @assert size(a) == dims
    @assert size(b) == dims
    @assert size(L) == dims
    @assert size(σ) == dims
    @assert size(η) == dims
end


## Various functions/methods for stiffness
"""
    shear_traction(::Type{<:PlaneFault}, u, λ, μ, dip)

Calculate the shear traction on the fault plane w.r.t. fault types.

## Arguments
- `u::AbstractArray{<:Number, 1}`: the output from [`dc3d_okada`](@ref)
- `λ::Number`: Lamé's first parameter
- `μ::Number`: shear modulus
- `dip::Number`: plane dip angle

## Reference
- A good reference is at [Displacement & Strain & Stress](https://nnakata.oucreate.com/page/Teaching_files/GEOPHYS130/GEOPHYS130_notes_all.pdf).
"""
shear_traction(ftype::Type{<:AbstractFault}) = NaN

@inline @traitfn function shear_traction(::Type{FT}, u::AbstractVector, λ::T, μ::T, dip::T) where {FT<:PlaneFault, T<:Number; IsDisplacementOnYZPlane{FT}}
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180))
end

@inline @traitfn function shear_traction(::Type{FT}, u::AbstractVector, λ::T, μ::T, dip::T) where {FT<:PlaneFault, T<:Number; !IsDisplacementOnYZPlane{FT}}
    # Special case for `dip = 90.` against above, however, at x-y plane instead.
    μ * (u[5] + u[7])
end

"""
    stiffness_tensor(fa::PlaneFaultDomain, gd::BoundaryElementGrid, ep::HomogeneousElasticProperties)

Calculate the reduced stiffness tensor. For 2D fault, the final result will be dimensionally reduced to a 3D array
    due to the *translational* & *reflective* & *perodic* symmetry, such that the tensor contraction will be equivalent to convolution,
    hence we could use FFT for better performace.

## Note
- Faults are originated from surface and extends downwards, thus `dep = 0`
"""
function stiffness_tensor(fa::PlaneFaultDomain{ftype, 1, T}, gd::BoundaryElementGrid{1}, ep::HomogeneousElasticProperties; kwargs...) where {ftype<:PlaneFault, T<:Number}
    st = SharedArray{T}(gd.nξ, gd.nξ)
    stiffness_tensor!(st, fa, gd, ep; kwargs...)
    return sdata(st)
end

function stiffness_tensor(fa::PlaneFaultDomain{ftype, 2, T}, gd::BoundaryElementGrid{2}, ep::HomogeneousElasticProperties; fourier_domain=true, kwargs...) where {ftype<:PlaneFault, T<:Number}
    st = SharedArray{T}(gd.nx, gd.nξ, gd.nξ)
    stiffness_tensor!(st, fa, gd, ep; kwargs...)

    function __convert_to_fourier_domain__()
        FFTW.set_num_threads(FFT_CONFIGS["FFT_NUM_THREADS"])
        x1 = zeros(T, 2 * gd.nx - 1)
        p1 = plan_rfft(x1, flags=FFT_CONFIGS["FFT_FLAG"])
        st_dft = Array{Complex{T}}(undef, gd.nx, gd.nξ, gd.nξ)
        @inbounds for l = 1: gd.nξ, j = 1: gd.nξ
            # The most tricky part to ensure correct FFT, see `dτ_dt!` for 2D case as well.
            # Ref -> (https://github.com/JuliaMatrices/ToeplitzMatrices.jl/blob/cbe29c344be8363f33eb17090121f8cff600b72e/src/ToeplitzMatrices.jl#L627)
            st_dft[:,j,l] .= p1 * [st[:,j,l]; reverse(st[2:end,j,l])]
        end
        return st_dft
    end

    fourier_domain ? __convert_to_fourier_domain__() : sdata(st)
end

function stiffness_chunk!(st::SharedArray{T, 2}, fa::PlaneFaultDomain{ftype, 1, T}, gd::BoundaryElementGrid{1}, ep::HomogeneousElasticProperties,
    subs::AbstractArray; kwargs...) where {ftype, T}

    unit_disl = applied_unit_dislocation(ftype)
    @inbounds @simd for sub in subs
        i, j = sub[1], sub[2]
        _u = dc3d_okada(zero(T), gd.y[i], gd.z[i], ep.α, zero(T), fa.dip, gd.ax, gd.aξ[j], unit_disl)
        st[i,j] = shear_traction(ftype, _u, ep.λ, ep.μ, fa.dip)
    end
end

function stiffness_chunk!(st::SharedArray{T, 3}, fa::PlaneFaultDomain{ftype, 2, T}, gd::BoundaryElementGrid{2}, ep::HomogeneousElasticProperties,
    subs::AbstractArray; nrept::Integer=2) where {ftype, T}

    unit_disl = applied_unit_dislocation(ftype)
    lrept = (gd.buffer_ratio + one(T)) * fa[:x]
    _u = Vector{T}(undef, 12)
    for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        # As stated in `BEMGrid_2D`, `depth` here are fixed at 0
        stiffness_periodic_boundary_condition!(_u, gd.x[i], gd.y[j], gd.z[j], ep.α, zero(T), fa.dip, gd.ax[1], gd.aξ[l], unit_disl, nrept, lrept)
        st[i,j,l] = shear_traction(ftype, _u, ep.λ, ep.μ, fa.dip)
    end
end

function stiffness_tensor!(st::SharedArray, fa::PlaneFaultDomain, gd::BoundaryElementGrid, ep::HomogeneousElasticProperties; kwargs...)
    @sync begin
        for p in procs(st)
            @async remotecall_wait(stiffness_shared_chunk!, p, st, fa, gd, ep; kwargs...)
        end
    end
end

function stiffness_shared_chunk!(st::SharedArray, fa::PlaneFaultDomain, gd::BoundaryElementGrid, ep::HomogeneousElasticProperties; kwargs...)
    i2s = CartesianIndices(st)
    inds = localindices(st)
    subs = i2s[inds]
    stiffness_chunk!(st, fa, gd, ep, subs; kwargs...)
end

"""
Periodic boundary condition for 2D faults.

## Arguments
- `u::AbstractVector`: In-place output which is a 12-elements vector (exactly the output of [`dc3d_okada`](@ref)). No assertion here imposed.
- same as [`dc3d_okada`](@ref), see [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) for details.
- `nrept::Integer`: (half) number of repetition, as denoted by `-npret: nrept`
- `lrept::Number`: length of repetition interval, see *Note* below

## Note
- The buffer block length is (`buffer_ratio` - 1) multipled by along-strike length.
"""
function stiffness_periodic_boundary_condition!(u::AbstractVector{T}, x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    nrept::Integer, lrept::T) where {T<:Number}
    fill!(u, zero(T))
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i * lrept, aw, disl)
    end
end

"""
For noraml fault, it should of course be [0., -1., 0.]. However, in term of force balance, it is quivalent to thrust fault
if dip angle are constrained within [0, π/2] in fact.

The unit of *unit dislocation* below is the same of `v * t` at set by user so to avoid normalization step.
"""
applied_unit_dislocation(::Type{NormalFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{ThrustFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{StrikeSlipFault}) = [1., 0., 0.]


## Core functions/methods for boundary element methods
"Temporal variable in solving ODEs aimed to avoid allocation overheads."
abstract type ODEStateVariable{dim} end

struct ODEState_1D{T<:AbstractVecOrMat{<:Real}} <: ODEStateVariable{1}
    dτ_dt::T
    relv::T
end

struct ODEState_2D{T<:AbstractArray{<:Real}, U<:AbstractArray{<:Complex}, P<:Plan} <: ODEStateVariable{2}
    dτ_dt::T
    relv::T
    dτ_dt_dft::U
    relv_dft::U
    dτ_dt_buffer::T
    pf::P
end

create_tmp_var(nξ::Integer; T1=Float64) = ODEState_1D([Vector{T1}(undef, nξ) for _ in 1: 2]...)

function create_tmp_var(nx::I, nξ::I; T1=Float64) where {I <: Integer}
    FFTW.set_num_threads(FFT_CONFIGS["FFT_NUM_THREADS"])
    x1 = Matrix{T1}(undef, 2 * nx - 1, nξ)
    p1 = plan_rfft(x1, 1, flags=FFT_CONFIGS["FFT_FLAG"])

    tvar = ODEState_2D(
        Matrix{T1}(undef, nx, nξ),
        zeros(T1, 2nx-1, nξ), # for relative velocity
        [Matrix{Complex{T1}}(undef, nx, nξ) for _ in 1: 2]...,
        Matrix{T1}(undef, 2nx-1, nξ),
        p1)
end

create_tmp_var(gd::BoundaryElementGrid{1}) = create_tmp_var(gd.nξ; T1=typeof(gd.Δξ))
create_tmp_var(gd::BoundaryElementGrid{2}) = create_tmp_var(gd.nx, gd.nξ; T1=typeof(gd.Δx))

function derivations!(du, u, mp::PlaneMaterialProperties{dim}, tvar::ODEStateVariable{dim}, se::StateEvolutionLaw, fform::FrictionLawForm) where {dim}
    _ndim = dim + 1
    v = selectdim(u, _ndim, 1)
    θ = selectdim(u, _ndim, 2)
    dv = selectdim(du, _ndim, 1)
    dθ = selectdim(du, _ndim, 2)

    variable_trim!(v, θ)
    dτ_dt!(mp, tvar, v)
    dv_dθ_dt!(fform, dv, dθ, v, θ, mp, tvar, se)
end

@inline function variable_trim!(v::T, θ::T) where {T<:AbstractVecOrMat{<:Float64}}
    clamp!(θ, 0., Inf)
end

@inline function dτ_dt!(mp::PlaneMaterialProperties{1}, tvar::ODEStateVariable{1}, v::AbstractArray)
    @fastmath @inbounds @simd for i = 1: prod(mp.dims)
        tvar.relv[i] = mp.vpl - v[i]
    end
    mul!(tvar.dτ_dt, mp.k, tvar.relv)
end

@inline function dτ_dt!(mp::PlaneMaterialProperties{2}, tvar::ODEStateVariable{2}, v::AbstractArray{T}) where {T<:Number}
    @fastmath @inbounds for j = 1: mp.dims[2]
        @simd for i = 1: mp.dims[1]
            tvar.relv[i,j] = mp.vpl - v[i,j]
        end
    end
    mul!(tvar.relv_dft, tvar.pf, tvar.relv)
    fill!(tvar.dτ_dt_dft, zero(T))

    @fastmath @inbounds for l = 1: mp.dims[2]
        for j = 1: mp.dims[2]
            @simd for i = 1: mp.dims[1]
                tvar.dτ_dt_dft[i,j] += mp.k[i,j,l] * tvar.relv_dft[i,l]
            end
        end
    end

    ldiv!(tvar.dτ_dt_buffer, tvar.pf, tvar.dτ_dt_dft)

    @fastmath @inbounds for j = 1: mp.dims[2]
        @simd for i = 1: mp.dims[1]
            tvar.dτ_dt[i,j] = tvar.dτ_dt_buffer[i,j]
        end
    end
end

@inline function dv_dθ_dt!(::CForm, dv::T, dθ::T, v::T, θ::T, mp::PlaneMaterialProperties, tvar::ODEStateVariable, se::StateEvolutionLaw) where {T<:AbstractVecOrMat}
    @fastmath @inbounds @simd for i = 1: prod(mp.dims)
        dμ_dθ = mp.σ[i] * mp.b[i] / θ[i]
        dμ_dv = mp.σ[i] * mp.a[i] / v[i]
        dθ[i] = dθ_dt(se, v[i], θ[i], mp.L[i])
        dv[i] = dv_dt(tvar.dτ_dt[i], dμ_dv, dμ_dθ, dθ[i], mp.η[i])
    end
end

@inline function dv_dθ_dt!(::RForm, dv::T, dθ::T, v::T, θ::T, mp::PlaneMaterialProperties, tvar::ODEStateVariable, se::StateEvolutionLaw) where {T<:AbstractVecOrMat}
    @fastmath @inbounds @simd for i = 1: prod(mp.dims)
        ψ1 = exp((mp.f0 + mp.b[i] * log(mp.v0 * θ[i] / mp.L[i])) / mp.a[i]) / 2mp.v0
        ψ2 = mp.σ[i] * ψ1 / hypot(1, v[i] * ψ1)
        dμ_dv = mp.a[i] * ψ2
        dμ_dθ = mp.b[i] / θ[i] * v[i] * ψ2
        dθ[i] = dθ_dt(se, v[i], θ[i], mp.L[i])
        dv[i] = dv_dt(tvar.dτ_dt[i], dμ_dv, dμ_dθ, dθ[i], mp.η[i])
    end
end

"""
    EarthquakeCycleProblem(p::PlaneMaterialProperties, u0, tspan; se=DieterichStateLaw(), fform=CForm())

Return an `ODEProblem` that encapsulate all the parameters and functions required for simulation. For the entailing usage, please refer [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)

## Arguments
- `gd::BoundaryElementGrid`: grids for fault domain.
- `p::PlaneMaterialProperties`: material profile.
- `u0::AbstractArray`: initial condition, should be organized such that the first of last dim is velocity while the 2nd of last dim is state.
- `tspan::NTuple`: time interval to be simulated.
- `se::StateEvolutionLaw`: state evolution law to be applied.
- `fform::FrictionLawForm`: forms of frictional law to be applied.
"""
function EarthquakeCycleProblem(gd::BoundaryElementGrid, p::PlaneMaterialProperties, u0::AbstractArray, tspan::NTuple; se=DieterichStateLaw(), fform=CForm()) where {dim}
    (fform == RForm() && minimum(p.η) ≈ 0.0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    tvar = create_tmp_var(gd)
    f! = (du, u, p, t) -> derivations!(du, u, p, tvar, se, fform)
    ODEProblem(f!, u0, tspan, p)
end

function friction(mp::PlaneMaterialProperties{N}, v::AbstractVecOrMat, θ::AbstractVecOrMat; fform=CForm()::FrictionLawForm) where {N}
    friction.(Ref(fform), v, θ, mp.L, mp.a, mp.b, mp.f0, mp.v0)
end

function friction(mp::PlaneMaterialProperties{N}, vθ::AbstractArray; kwargs...) where {N}
    ndim = ndims(vθ)
    friction(mp, selectdim(vθ, ndim, 1), selectdim(vθ, ndim, 2); kwargs...)
end
