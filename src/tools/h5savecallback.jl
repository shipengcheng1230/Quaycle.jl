export @h5savecallback

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
    V1<:AbstractVector, V2<:Tuple, V3<:Tuple, UT<:AbstractUnitRange}
    file::S
    buffer::D
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

function create_h5buffer(file::AbstractString, u::ArrayPartition, nstep::Integer, tstop::Real, ustrs, tstr)
    @assert tstr ∉ ustrs "Duplicate name of $(tstr) in $(ustrs)."
    @assert length(ustrs) == length(u.x) "Unmatched length between array partitions and names."
    buffer = h5savebufferzone(u, nstep, ustrs)
    buffer[tstr] = Vector{eltype(u)}(undef, nstep)
    count, total = 1, 0
    uiter = Base.OneTo(length(u.x))
    ushapes = map(size, u.x)
    f = x -> map(Base.Slice, axes(x))
    idxs = map(f, u.x)
    du = similar(u)
    h5open(file, "w") do f
        d_create(f, tstr, datatype(typeof(tstop)), ((nstep,), (-1,)), "chunk", (nstep,))
        for i ∈ uiter
            accusize = (ushapes[i]..., nstep)
            d_create(f, ustrs[i], datatype(eltype(u.x[i])), (accusize, (ushapes[i]..., -1,)), "chunk", accusize)
        end
    end
    return H5SaveBuffer(file, buffer, nstep, count, total, uiter, ustrs, ushapes, idxs, du, tstop, tstr)
end

h5savebufferzone(u::AbstractArray, nstep::Integer) = Array{eltype(u)}(undef, size(u)..., nstep)
h5savebufferzone(u::ArrayPartition, nstep, names) = h5savebufferzone(u.x, nstep, names)
h5savebufferzone(u::Tuple, nstep, names) = Dict(names[i] => h5savebufferzone(u[i], nstep) for i in 1: length(u))

function h5savebuffercbkernel(u, t, integrator, b::H5SaveBuffer, getu::Function)
    ptrs = getu(u, t, integrator, b)
    _trigger_copy(b, ptrs, t)
    (t == b.tstop || b.count > b.nstep) && _trigger_save(b, ptrs, t)
end

function _trigger_copy(b::H5SaveBuffer, ptrs, t)
    b.buffer[b.tstr][b.count[1]] = t
    for i ∈ b.uiter
        b.buffer[b.ustrs[i]][b.idxs[i]..., b.count] .= ptrs[i]
    end
    b.count += 1
end

function _trigger_save(b::H5SaveBuffer, ptrs, t)
    h5open(b.file, "r+") do f
        ht = d_open(f, b.tstr)
        set_dims!(ht, (b.total + b.count - 1,))
        ht[b.total+1: b.total+b.count-1] = b.count > b.nstep ? b.buffer[b.tstr] : selectdim(b.buffer[b.tstr], 1, 1: b.count-1)
        for i ∈ b.uiter
            hd = d_open(f, b.ustrs[i])
            set_dims!(hd, (b.ushapes[i]..., b.total + b.count - 1))
            hd[b.idxs[i]..., b.total+1: b.total+b.count-1] = b.count > b.nstep ? b.buffer[b.ustrs[i]] :
                view(b.buffer[b.ustrs[i]], b.idxs[i]..., 1: b.count-1)
        end
    end
    b.total += b.count - 1
    b.count = 1
end

function GET_VSSAR(u, t, integrator, buffer::H5SaveBuffer)
    get_du!(buffer.du, integrator)
    return (u.x[1], u.x[2], buffer.du.x[3])
end
