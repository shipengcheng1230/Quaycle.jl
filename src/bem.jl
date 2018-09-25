using SimpleTraits
using Parameters
using FFTW
using FFTW: Plan

using Distributed
using Base.Threads
using LinearAlgebra
using SharedArrays
using FileIO
using JLD2

include(joinpath(@__DIR__, "dc3d.jl"))

export discretize, properties, stiffness_tensor
export EarthquakeCycleProblem
export dc3d_okada

@traitdef IsOnYZPlane{faulttype}
@traitimpl IsOnYZPlane{faulttype} <- isonyz(faulttype)

isonyz(::Type{NormalFault}) = true
isonyz(::Type{ThrustFault}) = true
isonyz(::Type{StrikeSlipFault}) = false

abstract type AbstractDifferenceGrids{dim, isuniform, T} end
abstract type BoundaryElementGrid{dim, isuniform, T} <: AbstractDifferenceGrids{dim, isuniform, T} end

struct BEMGrid_1D{U1<:AbstractVector, U2<:AbstractMatrix, T<:Number, I<:Integer} <: BoundaryElementGrid{1, true, T}
    ξ::U1 # along-downdip
    Δξ::T
    nξ::I
    aξ::U2 # cell boundary along-downdip (width)
    y::U1 # y-component of `ξ`
    z::U1 # z-component of `ξ`
    ax::U1 # (psudo) cell boundary along-strike (length)
end

struct BEMGrid_2D{U1<:AbstractVector, U2<:AbstractMatrix, T<:Number, I<:Integer} <: BoundaryElementGrid{2, true, T}
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
    buffer::T # buffer on the edge of along-strike
end

"""
    discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ::T; ax_ratio=50*one(T)::T)

Generate the grid for given 1D fault domain.

## Arguments
- `Δξ`: grid space along-downdip
- `ax_ratio::Number`: ration of along-strike length agsinst along-downdip length for mimicing an extended
    2d (x & ξ) fault represented by 1d (ξ) domain. Default `ax_ratio=12.5` is more than enough for producing consistent results.
    Exceedingly large `ax_ratio` could potentially cause inconsistency.
"""
function discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ::T; ax_ratio=12.5) where {ftype, T}
    ξ, nξ, aξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    ax = fa[:ξ] * ax_ratio .* [-one(T), one(T)]
    BEMGrid_1D(ξ, Δξ, nξ, aξ, ξ.*cospi(fa.dip/180), ξ.*sinpi(fa.dip/180), ax)
end

"""
    discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ; buffer=:auto) where {ftype <: PlaneFault}

Generate the grid for given 2D fault domain.

## Arguments
- `Δx, Δξ`: grid space along-strike and along-downdip respectively
- `buffer::Union{T, Symbol}`: length of buffer size for introducing *zero-dislocation* area at along-strike edges of defined fault domain.
"""
function discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ; buffer=:auto) where {ftype <: PlaneFault}
    ξ, nξ, aξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    x, nx, ax = _divide_segment(Val(:symmatzero), fa[:x], Δx)
    if buffer == :auto
        buffersize = ftype == StrikeSlipFault ? fa[:x] : zero(typeof(fa[:x]))
    elseif typeof(buffer) <: Number
        buffersize = promote(buffer, fa[:x])[1]
    else
        error("Received buffer: $buffer. Valid `buffer` should be `:auto` or a number.")
        return nothing
    end
    BEMGrid_2D(x, ξ, Δx, Δξ, nx, nξ, ax, aξ, ξ.*cospi(fa.dip/180), ξ.*sinpi(fa.dip/180), buffersize)
end

discretize(fa::PlaneFaultDomain{ftype, 1}; nξ=500, kwargs...) where {ftype} = discretize(fa, fa[:ξ]/nξ; kwargs...)
discretize(fa::PlaneFaultDomain{ftype, 2}; nx=64, nξ=32, buffer=:auto) where  {ftype} = discretize(fa, fa[:x]/nx, fa[:ξ]/nξ; buffer=buffer)

function _divide_segment(::Val{:halfspace}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(zero(T), stop=-x+Δx, step=-Δx)) .- Δx/2
    xs = cat(xi .- Δx/2, xi .+ Δx/2, dims=2)
    return xi, length(xi), xs
end

function _divide_segment(::Val{:symmatzero}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx))
    xs = cat(xi .- Δx/2, xi .+ Δx/2, dims=2)
    return xi, length(xi), xs
end

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

@with_kw struct PlaneMaterialProperties{N, U<:AbstractVecOrMat, P<:AbstractArray, T<:Number} <: AbstractMaterialProperties{N}
    nsize::NTuple{N} # broadcasting size of array below
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    k::P # stiffness tensor
    σ::U # effective normal stress
    η::U # radiation damping
    vpl::T # plate rate, unlike pure Rate-State Friction simulation, here is restrained to be constant
    f0::T # ref. frictional coeff
    v0::T # ref. velocity
end

"""
    properties(fa::PlaneFaultDomain, gd::BoundaryElementGrid{dim}; _kwargs...) where {dim}

Establishing a material-properties-profile given by the fault domain and grids. User must provide the
    necessary parameters in according to the grid size specified or just a scalar for broadcasting.

## Arguments that are needed:
- `a`: contrib from velocity.
- `b`: contrib from state.
- `L`: critical distance.
- `σ`: effective normal stress.
- `vpl`: plate rate.
- `f0`: ref. frictional coeff.
- `v0`: ref. velocity.

## Arguments that are optional
- `k`: stiffness tensor. If `:auto`, it will automatically calculate by seeking `λ` and `μ` otherwise should be a valid file path to a `*.jld2` or an `AbstractArray`.
- `η`: radiation damping. If `:auto`, it will automatically seek `μ` and `vs` and use ``μ / 2\\mathrm{Vs}``.
- `vs`: shear wave velocity.
- `λ`: Lamé's first parameter
- `μ`: shear modulus
"""
function properties(fa::PlaneFaultDomain, gd::BoundaryElementGrid{dim}; _kwargs...) where {dim}

    function get_k()
        _k = get(kwargs, :k, nothing)
        _tk = typeof(_k)
        if _tk <: AbstractString && isfile(_k) && endswith(_k, ".jld2")
            @info "Loading stiffness: $_k ..."
            @load _k k
        elseif _tk <: AbstractArray
            k = _k
        elseif _k == :auto || _k == nothing
            λ = args_get_expand(:λ, kwargs, (), false)
            μ = args_get_expand(:μ, kwargs, (), false)
            ep = HomogeneousElasticProperties(λ=μ, μ=μ)
            @info "Calculating stiffness tensor..."
            k = stiffness_tensor(fa, gd, ep)
        else
            error("""
                Invalid option: $_k, should be:
                (1) A valid file path to a `*.jld2` file
                (2) An `AbstractArray`
                (3) Symbol `:auto`.
            """)
        end
        return k
    end

    function get_η()
        _η = get(kwargs, :η, nothing)
        if _η == :auto || _η == nothing
            μ = args_get_expand(:μ, kwargs, (), false)
            vs = args_get_expand(:vs, kwargs, gsize)
            η = μ ./ 2vs
        else
            η = args_get_expand(:η, kwargs, gsize)
        end
    end

    gsize = dim == 1 ? (gd.nξ,) : (gd.nx, gd.nξ)
    kwargs = _kwargs.data
    a = args_get_expand(:a, kwargs, gsize)
    b = args_get_expand(:b, kwargs, gsize)
    L = args_get_expand(:L, kwargs, gsize)
    σ = args_get_expand(:σ, kwargs, gsize)

    vpl = args_get_expand(:vpl, kwargs, (), false)
    f0 = args_get_expand(:f0, kwargs, (), false)
    v0 = args_get_expand(:v0, kwargs, (), false)

    k = get_k()
    η = get_η()

    @info "Establishing material properties..."
    mp = PlaneMaterialProperties(nsize=gsize, a=a, b=b, L=L, k=k, σ=σ, η=η, vpl=vpl, f0=f0, v0=v0)
    return mp
end

function args_get_expand(sym::Symbol, kwargs::NamedTuple, gsize::NTuple, expand::Bool=true)
    x = get(kwargs, sym, nothing)
    x == nothing && error("`$sym` is not provided.")
    if expand
        if typeof(x) <: Number
            x = x .* ones(gsize...)
        elseif typeof(x) <: AbstractArray
            size(x) == gsize || error("Dim `$sym` does not match with given grid, $(size(x)) received, $gsize required.")
        else
            error("Illegal input type of `$sym`, get $x, required `Number` or `AbstractArray`.")
        end
    end
    return x
end

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

@inline @traitfn function shear_traction(::Type{FT}, u, λ, μ, dip) where {FT<:PlaneFault; IsOnYZPlane{FT}}
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180))
end

@inline @traitfn function shear_traction(::Type{FT}, u, λ, μ, dip) where {FT<:PlaneFault; !IsOnYZPlane{FT}}
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
function stiffness_tensor(fa::PlaneFaultDomain, gd::BoundaryElementGrid, ep::HomogeneousElasticProperties; kwargs...)
end

function stiffness_tensor(fa::PlaneFaultDomain{ftype, 1, T}, gd::BoundaryElementGrid{1, true}, ep::HomogeneousElasticProperties,
    ) where {ftype<:PlaneFault, T<:Number}

    udisl = applied_unit_dislocation(ftype)
    if nprocs() == 1
        u12 = zeros(T, 12)
        ST = zeros(T, gd.nξ, gd.nξ)
        for j = 1: gd.nξ
            @inbounds @simd for i = 1: gd.nξ
                u12 .= dc3d_okada(zero(T), gd.y[i], gd.z[i], ep.α, zero(T), fa.dip, gd.ax, gd.aξ[j,:], udisl)
                ST[i,j] = shear_traction(ftype, u12, ep.λ, ep.μ, fa.dip)
            end
        end
        return ST
    else
        return stiffness_tensor_parfor(fa, gd, ep, udisl)
    end
end

function stiffness_tensor_parfor(fa::PlaneFaultDomain{ftype, 1, T}, gd::BoundaryElementGrid{1, true}, ep::HomogeneousElasticProperties,
    _udisl::AbstractVector) where {ftype<:PlaneFault, T<:Number}

    ST = SharedArray{T}(gd.nξ, gd.nξ)
    @sync begin
        for p in procs(ST)
            @async remotecall_wait(()->(fa, gd, ep, _udisl), p)
        end
    end

    @sync @distributed for j = 1: gd.nξ
        @inbounds @simd for i = 1: gd.nξ
            _u = dc3d_okada(zero(T), gd.y[i], gd.z[i], ep.α, zero(T), fa.dip, gd.ax, gd.aξ[j,:], _udisl)
            ST[i,j] = shear_traction(ftype, _u, ep.λ, ep.μ, fa.dip)
        end
    end

    clear!([:fa, :gd, :ep, :_udisl, :_u], procs(ST))
    return sdata(ST)
end

function stiffness_tensor(fa::PlaneFaultDomain{ftype, 2, T}, gd::BoundaryElementGrid{2, true}, ep::HomogeneousElasticProperties;
    nrept::Integer=2) where {ftype<:PlaneFault, T<:Number}

    udisl = applied_unit_dislocation(ftype)
    ST = SharedArray{T}(gd.nx, gd.nξ, gd.nξ)

    @sync begin
        for p in procs(ST)
            @async remotecall_wait(stiffness_shared_chunk!, p, ST, fa, gd, ep, udisl, nrept)
        end
    end

    ST = sdata(ST)
    FFTW.set_num_threads(nthreads())
    x1 = zeros(T, gd.nx)
    p1 = plan_rfft(x1, flags=FFTW.MEASURE)
    ST_DFT = Array{Complex{T}}(undef, floor(Int, gd.nx/2)+1, gd.nξ, gd.nξ)
    @inbounds for l = 1: gd.nξ, j = 1: gd.nξ
        ST_DFT[:,j,l] .= p1 * ST[:,j,l]
    end
    return ST_DFT
end

function stiffness_chunk!(ST::SharedArray, fa::PlaneFaultDomain{ftype, 2, T}, gd::BoundaryElementGrid{2, true}, ep::HomogeneousElasticProperties,
    _udisl, nrept, subs) where {ftype, T}
    _u = zeros(T, 12)
    for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        _u .= stiffness_periodic_boundary_condition(gd.x[i], gd.y[j], gd.z[j], ep.α, zero(T), fa.dip, gd.ax[1,:], gd.aξ[l,:], _udisl, nrept, fa[:x]+gd.buffer)
        ST[i,j,l] = shear_traction(ftype, _u, ep.λ, ep.μ, fa.dip)
    end
end

function stiffness_shared_chunk!(ST::SharedArray, fa::PlaneFaultDomain{ftype, 2, T}, gd::BoundaryElementGrid{2, true}, ep::HomogeneousElasticProperties,
    _udisl, nrept) where {ftype, T}
    i2s = CartesianIndices(ST)
    inds = localindices(ST)
    subs = i2s[inds]
    stiffness_chunk!(ST, fa, gd, ep, _udisl, nrept, subs)
end

"""
Periodic boundary condition for 2D faults.

## Arguments
- same as `dc3d_okada`, see [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) for details.
- `ax::AbstractVector`: along-strike fault length
- `nrept::Integer`: (half) number of repetition, as denoted by `-npret: nrept`
- `lrept::Number`: length of repetition interval, see *Note* below

## Note
- The buffer block is evenly distributed on the two along-strike edges, each of which contains half of that.
"""
function stiffness_periodic_boundary_condition(x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    nrept::Integer, lrept::T) where {T<:Number}
    u = zeros(T, 12)
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i * lrept, aw, disl)
    end
    return u
end

"""
For noraml fault, it should of course be [0., -1., 0.]. However, in term of force balance, it is quivalent to thrust fault
if dip angle are constrained within [0, π/2] in fact.

The unit of *unit dislocation* below is the same of `v * t` at set by user so to avoid normalization step.
"""
applied_unit_dislocation(::Type{NormalFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{ThrustFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{StrikeSlipFault}) = [1., 0., 0.]

"Temporal variable in solving ODEs aimed to avoid allocation overheads."
struct TmpVariable{T<:AbstractVecOrMat{<:Real}, U<:Union{AbstractVecOrMat{<:Complex}, Nothing}, P1<:Union{Plan, Nothing}, P2<:Union{Plan, Nothing}}
    dτ_dt::T
    dμ_dθ::T
    dμ_dv::T
    ψ1::T
    ψ2::T
    relv::T
    relv_dft::U
    dτ_dt_dft::U
    plan_1d::P1
    plan_2d::P2
end

function create_tmp_var(nsize::NTuple{N}; T1=Float64) where {N}
    if N == 1
        tvar = TmpVariable([Vector{T1}(undef, nsize[1]) for _ in 1: 6]..., fill(nothing, 4)...)
    elseif N == 2
        FFTW.set_num_threads(nthreads())
        x1 = Vector{T1}(undef, nsize[1])
        x2 = Matrix{T1}(undef, nsize...)
        p1 = plan_rfft(x1, flags=FFTW.MEASURE)
        p2 = plan_rfft(x2, 1, flags=FFTW.MEASURE)
        rfft_size = floor(Int, nsize[1]/2) + 1
        tvar = TmpVariable(
            [Matrix{T1}(undef, nsize...) for _ in 1: 6]...,
            [Matrix{Complex{T1}}(undef, rfft_size, nsize[2]) for _ in 1: 2]...,
            p1, p2)
    else
        error("Dim: $N received. Valid `N` is `1` or `2`.")
    end
    return tvar
end

function derivations!(du, u, mp::PlaneMaterialProperties{dim}, tvar::TmpVariable, se::StateEvolutionLaw, fform::FrictionLawForm) where {dim}
    _ndim = dim + 1
    v = selectdim(u, _ndim, 1)
    θ = selectdim(u, _ndim, 2)
    dv = selectdim(du, _ndim, 1)
    dθ = selectdim(du, _ndim, 2)
    clamp!(θ, 0.0, Inf)

    @. tvar.relv = mp.vpl - v
    dθ_dt!(dθ, se, v, θ, mp.L)
    dτ_dt!(mp, tvar)
    dv_dθ_dt!(fform, dv, dθ, v, θ, mp, tvar)
end

@inline function dτ_dt!(mp::PlaneMaterialProperties{1}, tvar::TmpVariable)
    mul!(tvar.dτ_dt, mp.k, tvar.relv)
end

@inline function dτ_dt!(mp::PlaneMaterialProperties{2}, tvar::TmpVariable)
    tvar.relv_dft .= tvar.plan_2d * tvar.relv
    fill!(tvar.dτ_dt_dft, 0.)
    nd = size(tvar.dτ_dt, 2)
    @inbounds for j = 1: nd
        @simd for l = 1: nd
            tvar.dτ_dt_dft[:,j] .+= mp.k[:,j,l] .* tvar.relv_dft[:,l]
        end
    end
    ldiv!(tvar.dτ_dt, tvar.plan_2d, tvar.dτ_dt_dft)
end

dθ_dt!(dθ::T, se::StateEvolutionLaw, v::T, θ::T, L::U) where {T<:AbstractVecOrMat, U<:AbstractVecOrMat} = dθ .= dθ_dt.(Ref(se), v, θ, L)

@inline function dv_dθ_dt!(::CForm, dv::T, dθ::T, v::T, θ::T, mp::PlaneMaterialProperties, tvar::TmpVariable) where {T<:AbstractVecOrMat}
    @. tvar.dμ_dθ = mp.σ * mp.b / θ
    @. tvar.dμ_dv = mp.σ * mp.a / v
    @. dv = dv_dt(tvar.dτ_dt, tvar.dμ_dv, tvar.dμ_dθ, dθ, mp.η)
end

@inline function dv_dθ_dt!(::RForm, dv::T, dθ::T, v::T, θ::T, mp::PlaneMaterialProperties, tvar::TmpVariable) where {T<:AbstractVecOrMat}
    @. tvar.ψ1 = exp((mp.f0 + mp.b * log(mp.v0 * θ / mp.L)) / mp.a) / 2mp.v0
    @. tvar.ψ2 = mp.σ * tvar.ψ1 / hypot(1, v * tvar.ψ1)
    @. tvar.dμ_dv = mp.a * tvar.ψ2
    @. tvar.dμ_dθ = mp.b / θ * v * tvar.ψ2
    @. dv = dv_dt(tvar.dτ_dt, tvar.dμ_dv, tvar.dμ_dθ, dθ, mp.η)
end

"""
    EarthquakeCycleProblem(p::PlaneMaterialProperties, u0, tspan; se=DieterichStateLaw(), fform=CForm())

Return an `ODEProblem` that encapsulate all the parameters and functions required for simulation. For the entailing usage, please refer [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)

## Arguments
- `p::PlaneMaterialProperties`: material profile.
- `u0::AbstractArray`: initial condition, should be organized such that the first of last dim is velocity while the 2nd of last dim is state.
- `tspan::NTuple`: time interval to be simulated.
- `se::StateEvolutionLaw`: state evolution law to be applied.
- `fform::FrictionLawForm`: forms of frictional law to be applied.
"""
function EarthquakeCycleProblem(p::PlaneMaterialProperties{dim}, u0::AbstractArray, tspan::NTuple; se=DieterichStateLaw(), fform=CForm()) where {dim}
    (fform == RForm() && minimum(p.η) ≈ 0.0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    tvar = create_tmp_var(p.nsize)
    f! = (du, u, p, t) -> derivations!(du, u, p, tvar, se, fform)
    ODEProblem(f!, u0, tspan, p)
end
