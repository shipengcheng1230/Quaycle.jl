abstract type AbstracGreensFunction end
abstract type AbstractAllocation{dim} end

"Obtain mapping from local linear index to cartesian index."
function get_subs(A::SharedArray)
    i2s = CartesianIndices(A)
    inds = localindices(A)
    subs = i2s[inds]
end

"A general pattern parallel SharedArray computation. See [here](https://docs.julialang.org/en/v1/manual/parallel-computing/#man-shared-arrays-1) for an example."
macro gen_shared_chunk_call(name::Symbol, isvector)
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
        if $(!isvector)
            function $(chunkstr)(st::SharedArray, args...; kwargs...)
                subs = get_subs(st)
                $(chunkstr)(st, subs, args...; kwargs...)
            end
        else
            function $(chunkstr)(st::NTuple{N, <:SharedArray}, args...; kwargs...) where N
                subs = get_subs(st[1]) # all SharedArray shall have the same dims
                $(chunkstr)(st, subs, args...; kwargs...)
            end
        end
    end)
end

## helper function
heaviside(x::T) where T = x ≤ zero(T) ? zero(T) : one(T)

## kernel function
const KERNELDIR = joinpath(@__DIR__, "gf")
foreach(x -> include(joinpath(KERNELDIR, x)), filter!(x -> endswith(x, ".jl"), readdir(KERNELDIR)))

## concrete greens function
include("gf_okada.jl")
