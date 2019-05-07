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
                    @async remotecall_wait($(chunkstr), WorkerPool(workers()), st, args...; kwargs...)
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

include("gf_okada.jl")
