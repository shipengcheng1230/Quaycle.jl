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
include("gf_okada.jl")
include("gf_sbarbot.jl")
include("gf_operator.jl")

## composite green's function
struct ViscoelasticCompositeGreensFunction{T1, T2, T3, T4}
    ee::T1 # elastic ⟷ elastic
    ev::T2 # elastic ⟷ inelastic
    ve::T3 # inelastic ⟷ elastic
    vv::T4 # inelastic ⟷ inelastic
end
