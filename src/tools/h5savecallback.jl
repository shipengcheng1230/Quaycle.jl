export @h5savecallback

"""
    @h5savecallback(filename, tend, nsteps, usize, T)

Construct a `FunctionCallingCallback` for incrementally stored output into HDF5 file.

## Arguments
- `filename`: file name to be stored
- `tend`: end time of simulation
- `nsteps`: after `nsteps` steps, a saving operation is performed otherwise caching them
- `T`: type of stored data
"""
macro h5savecallback(filename, tend, nsteps, usize, T)
    callback = gensym(:callback)
    nd = eval(:(length($(usize))))

    # Well, `esc` the whole expr might be a bad idea (https://discourse.julialang.org/t/macros-and-modules/5859/2)
    # but all the function called here are reexport from "HDF5", "DiffEqCallbacks" and "Base"
    # so no cross module mess
    # Besides, the whole `esc` returns crystal clear readable code by `@macroexpand` if you would like to probe.
    esc(quote
        let count = 1
            global $(callback)
            accu = Array{$T}(undef, $(usize.args...), $(nsteps))
            accusize = tuple($(usize.args...), $(nsteps))
            acct = Vector{$T}(undef, $(nsteps))
            acctsize = ($(nsteps),)
            total = 0
            h5open($(filename), "w") do fid
                d = d_create(fid, "u", $(T), (accusize, ntuple(_ -> -1, Val($(nd+1)))), "chunk", accusize)
                d = d_create(fid, "t", $(T), (acctsize, (-1,)), "chunk", acctsize)
            end

            function $(callback)(u, t, integrator)
                if t == $(tend)
                    rest = total % $(nsteps)
                    selectdim(accu, $(nd+1), count) .= u
                    selectdim(acct, 1, count) .= t
                    h5open($filename, "r+") do f
                        d = d_open(f, "u")
                        d[$((:(:) for _ in 1: nd)...), total-rest+1: total+1] = selectdim(accu, $(nd+1), 1: rest+1)
                        set_dims!(d, ($(usize.args...), total+1))

                        ht = d_open(f, "t")
                        ht[total-rest+1: total+1] = selectdim(acct, 1, 1: rest+1)
                        set_dims!(ht, (total+1,))
                    end
                elseif count > $(nsteps)
                    h5open($filename, "r+") do f
                        d = d_open(f, "u")
                        d[$((:(:) for _ in 1: nd)...), total-$(nsteps-1): total] = accu
                        set_dims!(d, ($(usize.args...), total+$(nsteps)))

                        ht = d_open(f, "t")
                        ht[total-$(nsteps-1): total] = acct
                        set_dims!(ht, (total+$(nsteps),))
                    end
                    selectdim(accu, $(nd+1), 1) .= u
                    selectdim(acct, 1, 1) .= t
                    count = 2
                    total += 1
                else
                    selectdim(accu, $(nd+1), count) .= u
                    selectdim(acct, 1, count) .= t
                    count += 1
                    total += 1
                end
            end

            FunctionCallingCallback($(callback); func_everystep=true)
        end
    end)
end
