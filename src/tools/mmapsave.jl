export mmapsave

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
