export mmapsave, @h5save, @setlastindex

"""
    DECallbackSaveToFile(iot::IOStream, iou::IOStream)

Construct a functional callback to write `ODESolution` (`t` & `u`) into file. The reason to separate `t` and `u` is
for more easily reshape `u` w.r.t grids specification. It right now falls on users' memory on what the type of solution is
for accurately retrieving results.

## Arguments
- `iot::IOStream`: stream pointing to solution of time
- `iou::IOStream`: stream pointing to solution of domain

**Note**
It is strongly not recommended to use "skipping" scheme (by defining `thrd` and `dts(a)` for each case) when solution is too
oscillated.
"""
function mmapsave(iot::IOStream, iou::IOStream)

    function __save(u, t, integrator)
        write(iot, t)
        write(iou, u)
    end
    FunctionCallingCallback(__save; func_everystep=true)
end

function mmapsave(iot::IOStream, iou::IOStream, thrd::Number; dts=0.1, dta=3.1536e6)

    function __save(u::AbstractArray{T, N}, t, integrator) where {T, N}
        maxv = maximum(selectdim(u, N, 1))
        maxvprev = maximum(selectdim(integrator.uprev, N, 1))
        if (maxv ≥ thrd && t - integrator.tprev ≥ dts) || (maxvprev < thrd && t - integrator.tprev ≥ dta)
            write(iot, t)
            write(iou, u)
        end
    end
    FunctionCallingCallback(__save; func_everystep=true)
end

macro setlastindex(dest, src, nd, dest_index)
    Expr(
        :(=),
        Expr(:ref, esc(dest), [:(:) for _ in 1: nd]..., esc(dest_index)),
        esc(src),
    )
end

macro setlastindex(dest, src, nd, dest_index, src_index)
    Expr(
        :(=),
        Expr(:ref, esc(dest), [:(:) for _ in 1: nd]..., esc(dest_index)),
        Expr(:call, :selectdim, esc(src), nd + 1, esc(src_index)),
    )
end

macro h5save(filename, tend, nsteps, usize, nd, T)
    callback = gensym(:callback)
    esc(quote
        let count = 1
            global $(callback)
            accu = Array{$T}(undef, $(usize)..., $(nsteps))
            accusize = tuple($((usize))..., $(nsteps))
            acct = Vector{$T}(undef, $(nsteps))
            acctsize = ($(nsteps),)
            total = 0
            h5open($(filename), "w") do fid
                d = d_create(fid, "u", $(T), (accusize, ntuple(_ -> -1, Val($nd+1))), "chunk", accusize)
                d = d_create(fid, "t", $(T), (acctsize, (-1,)), "chunk", acctsize)
            end

            function $callback(u, t, integrator)
                if t == $(tend)
                    rest = total % $(nsteps)
                    selectdim(accu, $nd+1, count) .= u
                    selectdim(acct, 1, count) .= t
                    h5open($filename, "r+") do f
                        d = d_open(f, "u")
                        @setlastindex(d, accu, $nd, total-rest+1: total+1, 1: rest+1)
                        set_dims!(d, ($(usize)..., total+1))

                        ht = d_open(f, "t")
                        @setlastindex(ht, acct, 0, total-rest+1: total+1, 1: rest+1)
                        set_dims!(ht, (total+1,))
                    end
                elseif count > $((nsteps))
                    h5open($filename, "r+") do f
                        d = d_open(f, "u")
                        @setlastindex(d, accu, $nd, total-$(nsteps)+1: total)
                        set_dims!(d, ($(usize)..., total+$(nsteps)))

                        ht = d_open(f, "t")
                        @setlastindex(ht, acct, 0, total-$(nsteps)+1: total)
                        set_dims!(ht, (total+$(nsteps),))
                    end
                    selectdim(accu, $nd+1, 1) .= u
                    selectdim(acct, 1, 1) .= t
                    count = 2
                    total += 1
                else
                    selectdim(accu, $nd+1, count) .= u
                    selectdim(acct, 1, count) .= t
                    count += 1
                    total += 1
                end
            end
        end
    end)
end
