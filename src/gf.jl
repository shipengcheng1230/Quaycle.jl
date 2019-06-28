export cat_greensfunc, compose_stress_greensfunc

## Static Green's Function
"Obtain mapping from local linear index to cartesian index."
function get_subs(A::SharedArray)
    i2s = CartesianIndices(A)
    inds = localindices(A)
    subs = i2s[inds]
end

"A general pattern parallel SharedArray computation. See [here](https://docs.julialang.org/en/v1/manual/parallel-computing/#man-shared-arrays-1) for an example."
macro gen_shared_chunk_call(name::Symbol)
    namestr = String(name) * "!"
    chunkstr = Symbol(namestr[1:end-1] * "_chunk!")
    namestr = Symbol(namestr)
    esc(quote
        function $(namestr)(st::S, args...; kwargs...) where S<:Union{SharedArray, NTuple{N, <:SharedArray}} where N
            pcs = typeof(st) <: SharedArray ? procs(st) : procs(st[1])
            @sync begin
                for p in pcs
                    @async remotecall_wait($(chunkstr), p, st, args...; kwargs...)
                end
            end
        end
        function $(chunkstr)(st::SharedArray, args...; kwargs...)
            subs = get_subs(st)
            $(chunkstr)(st, subs, args...; kwargs...)
        end
        function $(chunkstr)(st::NTuple{N, <:SharedArray}, args...; kwargs...) where N
            subs = get_subs(st[1]) # all SharedArray shall have the same dims
            $(chunkstr)(st, subs, args...; kwargs...)
        end
    end)
end

## helper function
heaviside(x::T) where T = x ≤ zero(T) ? zero(T) : one(T)
xlogy(x::T, y::T) where T = isapprox(x, zero(T)) ? zero(T) : x * log(y)
xlogy(x, y) = xlogy(promote(x, y)...)
omega(x::T) where T = heaviside(x + 1/2) - heaviside(x - 1/2)
S(x::T) where T = omega(x - 1/2)

## kernel function
const KERNELDIR = joinpath(@__DIR__, "gf")
foreach(x -> include(joinpath(KERNELDIR, x)), filter!(x -> endswith(x, ".jl"), readdir(KERNELDIR)))

## concrete greens function
export stress_greens_func
@gen_shared_chunk_call stress_greens_func

include("gf_dislocation.jl")
include("gf_strain.jl")
include("gf_operator.jl")

## composite green's function
struct ViscoelasticCompositeGreensFunction{T<:AbstractMatrix, U<:AbstractArray}
    ee::U # elastic ⟷ elastic
    ev::T # elastic ⟷ inelastic
    ve::T # inelastic ⟷ elastic
    vv::T # inelastic ⟷ inelastic
end

"""
    cat_greensfunc(ee::AbstractArray, ev::NTuple, ve::NTuple, vv::NTuple)

Concatenate tuple of matrix or tuple of tuple of matrix which arrange them in
    a way to update traction/stress rate using only one BLAS call. It does nothing
    to the elastic Green's function and specifically used for the outputs from
    [`stress_greens_func`](@ref) and [`stress_greens_func`](@ref).

## Arguments
- `ee`: traction Green's function within the elastic fault
- `ev`: stress Green's function from elastic fault to inelastic asthenosphere
- `ve`: traction Green's function inelastic asthenosphere to elastic fault
- `vv`: stress Green's function within inelastic asthenosphere
"""
function cat_greensfunc(ee::AbstractArray, ev::NTuple, ve::NTuple, vv::NTuple)
    σindex = Base.OneTo(length(ev))
    ϵindex = Base.OneTo(length(ve))
    m, n = size(ev[1]) # number of inelastic elements, number of fault patches
    ev′ = PseudoBlockArray{eltype(ev[1])}(undef, [m for _ in σindex], [n])
    ve′ = PseudoBlockArray{eltype(ve[1])}(undef, [n], [m for _ in ϵindex])
    vv′ = PseudoBlockArray{eltype(vv[1][1])}(undef, [m for _ in σindex], [m for _ in ϵindex])
    foreach(i -> setblock!(ev′, ev[i], i, 1), σindex)
    foreach(i -> setblock!(ve′, ve[i], 1, i), ϵindex)
    foreach(x -> setblock!(vv′, vv[x[2]][x[1]], x[1], x[2]), Iterators.product(σindex, ϵindex))
    ViscoelasticCompositeGreensFunction(ee, Array(ev′), Array(ve′), Array(vv′))
end

"""
    compose_stress_greensfunc(mf::OkadaMesh, me::SBarbotMeshEntity,
        λ::T, μ::T, ft::PlaneFault, comp::NTuple{N, <:Symbol}) where {T, N}

Shortcut function for computing all 4 Green's function for viscoelastic relaxation.
    Arguments stays the same as [`stress_greens_func`](@ref).
"""
function compose_stress_greensfunc(mf::OkadaMesh, me::SBarbotMeshEntity, λ::T, μ::T, ft::PlaneFault, comp::NTuple{N, <:Symbol}) where {T, N}
    ee = stress_greens_func(mf, λ, μ, ft)
    ev = stress_greens_func(mf, me, λ, μ, ft)
    ve = stress_greens_func(me, mf, λ, μ, ft, comp)
    vv = stress_greens_func(me, λ, μ, comp)
    return cat_greensfunc(ee, ev, ve, vv)
end
