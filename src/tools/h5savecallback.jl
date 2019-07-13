export @h5savecallback, wsolve, 𝐕𝚯, 𝐕𝚯𝚬, 𝐕𝚯𝚬′, 𝐕𝚯𝚬𝚺, 𝐕𝚬′

"""
    @h5savecallback(filename, tend, nsteps, usize, T)

Construct a `FunctionCallingCallback` for incrementally stored output into HDF5 file.
    This callback function only works for naive output arrays whose shape look like `A[..., :, :, :, ...]`.
    It is suggested to use this macro at top-level scope since it contains `eval`.

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

mutable struct H5SaveBuffer{
    S<:AbstractString, D<:AbstractDict, I<:Integer, A<:AbstractArray, R<:Real,
    V<:AbstractVector, V1<:AbstractVector, V2<:Tuple, V3<:Tuple, UT<:AbstractUnitRange}
    file::S
    ubuffer::D
    tbuffer::V
    nstep::I
    count::I
    total::I
    uiter::UT
    ustrs::V1
    ushapes::V2
    idxs::V3
    du::A
    tstop::R
    tstr::S
end

function create_h5buffer(file::AbstractString, ptrs::Tuple, du::AbstractArray, nstep::Integer, tstop::Real, ustrs, tstr)
    @assert tstr ∉ ustrs "Duplicate name of $(tstr) in $(ustrs)."
    @assert length(ustrs) == length(ptrs) "Unmatched length between solution components and names."
    ubuffer = h5savebufferzone(ptrs, nstep, ustrs)
    tbuffer = Vector{eltype(du)}(undef, nstep)
    count, total = 1, 0
    uiter = Base.OneTo(length(ptrs))
    ushapes = map(size, ptrs)
    f = x -> map(Base.Slice, axes(x))
    idxs = map(f, ptrs)
    h5open(file, "w") do f
        d_create(f, tstr, datatype(typeof(tstop)), ((nstep,), (-1,)), "chunk", (nstep,))
        for i ∈ uiter
            accusize = (ushapes[i]..., nstep)
            d_create(f, ustrs[i], datatype(eltype(du)), (accusize, (ushapes[i]..., -1,)), "chunk", accusize)
        end
    end
    return H5SaveBuffer(file, ubuffer, tbuffer, nstep, count, total, uiter, ustrs, ushapes, idxs, du, tstop, tstr)
end

h5savebufferzone(u::AbstractArray, nstep::Integer) = Array{eltype(u)}(undef, size(u)..., nstep)
h5savebufferzone(u::Tuple, nstep, names) = Dict(names[i] => h5savebufferzone(u[i], nstep) for i in 1: length(u))

function h5savebuffercbkernel(u, t, integrator, b::H5SaveBuffer, getu::Function)
    ptrs = getu(u, t, integrator)
    _trigger_copy(b, ptrs, t)
    (t == b.tstop || b.count > b.nstep) && _trigger_save(b, ptrs, t)
end

function _trigger_copy(b::H5SaveBuffer, ptrs, t)
    b.tbuffer[b.count[1]] = t
    for i ∈ b.uiter
        b.ubuffer[b.ustrs[i]][b.idxs[i]..., b.count] .= ptrs[i]
    end
    b.count += 1
end

function _trigger_save(b::H5SaveBuffer, ptrs, t)
    h5open(b.file, "r+") do f
        ht = d_open(f, b.tstr)
        set_dims!(ht, (b.total + b.count - 1,))
        ht[b.total+1: b.total+b.count-1] = ifelse(b.count > b.nstep, b.tbuffer, selectdim(b.tbuffer, 1, 1: b.count-1))
        for i ∈ b.uiter
            hd = d_open(f, b.ustrs[i])
            set_dims!(hd, (b.ushapes[i]..., b.total + b.count - 1))
            hd[b.idxs[i]..., b.total+1: b.total+b.count-1] = ifelse(b.count > b.nstep, b.ubuffer[b.ustrs[i]],
                view(b.ubuffer[b.ustrs[i]], b.idxs[i]..., 1: b.count-1))
        end
    end
    b.total += b.count - 1
    b.count = 1
end

# https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/785
"Retrieve **velocity**, **state** and **strain rate**."
𝐕𝚯𝚬′(u::ArrayPartition, t, integrator) = (u.x[1], u.x[2], integrator(integrator.t, Val{1}).x[3])
"Retrieve **velocity**, **state**, **strain** and **stress**."
𝐕𝚯𝚬𝚺(u::ArrayPartition, args...) = (u.x[1], u.x[2], u.x[3], u.x[4])
"Retrieve **velocity**, **state** and **strain**."
𝐕𝚯𝚬(u::ArrayPartition, args...) = (u.x[1], u.x[2], u.x[3])
"Retrieve **velocity** and **state**."
𝐕𝚯(u::ArrayPartition, args...) = (u.x[1], u.x[2])
"Retrieve **velocity** and **strain rate**"
𝐕𝚬′(u::ArrayPartition, t, integrator) = (u.x[1], integrator(integrator.t, Val{1}).x[3])

"""
    wsolve(prob::ODEProblem, alg::OrdinaryDiffEqAlgorithm,
        file, nstep, getu, ustrs, tstr; kwargs...)

Write the solution to HDF5 file while solving the ODE. The interface
    is exactly the same as
    [`solve` an `ODEProblem`](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html)
    except a few more about the saving procedure. Notice, it will set
    `save_everystep=false` so to avoid memory blow up. The return code
    will be written as an attribute in `tstr` data group.

## Extra Arguments
- `file::AbstractString`: name of file to be saved
- `nstep::Integer`: number of steps after which a saving operation will be performed
- `getu::Function`: function handler to extract desired solution for saving
- `ustr::AbstractVector`: list of names to be assigned for each components, whose
    length must equal the length of `getu` output
- `tstr::AbstractString`: name of time data
"""
function wsolve(prob::ODEProblem, alg::OrdinaryDiffEqAlgorithm, file, nstep, getu, ustrs, tstr; kwargs...)
    integrator = init(prob, alg)
    du = similar(prob.u0)
    ptrs = getu(prob.u0, prob.tspan[1], integrator)
    bf = create_h5buffer(file, ptrs, du, nstep, prob.tspan[2], ustrs, tstr)
    cb = (u, t, integrator) -> h5savebuffercbkernel(u, t, integrator, bf, getu)
    fcb = FunctionCallingCallback(cb)
    sol = solve(prob, alg; save_everystep=false, callback=fcb, kwargs...)
    _trigger_save(bf, ptrs, 0.0) # in case `solve` terminates earlier
    h5writeattr(file, tstr, Dict("retcode" => string(sol.retcode)))
    return sol
end
