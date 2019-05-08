abstract type AbstracGreensFunction end
abstract type AbstractAllocation{dim} end

macro gen_shared_chunk_call(name::Symbol)
    namestr = String(name) * "!"
    chunkstr = Symbol(namestr[1:end-1] * "_chunk!")
    namestr = Symbol(namestr)
    esc(quote
        function $(namestr)(st::SharedArray, args...; kwargs...) where T
            @sync begin
                for p in procs(st)
                    @async remotecall_wait($(chunkstr), p, st, args...; kwargs...)
                end
            end
        end

        function $(chunkstr)(st::SharedArray, args...; kwargs...) where T
            i2s = CartesianIndices(st)
            inds = localindices(st)
            subs = i2s[inds]
            $(chunkstr)(st, subs, args...; kwargs...)
        end
    end)
end

## helper function
heaviside(x::T) where T = x ≤ zero(T) ? zero(T) : one(T)

## kernel function
const KERNELDIR = joinpath(@__DIR__, "gfkernels")

for f in filter!(x -> endswith(x, ".jl"), readdir(KERNELDIR))
    include(abspath(joinpath(KERNELDIR, f)))
end

## concrete greens function
include("gf_okada.jl")
