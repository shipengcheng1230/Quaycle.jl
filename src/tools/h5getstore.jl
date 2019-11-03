export @getprop, @store, store

"""
    @getprop filename

Read property stored in HDF5.

## Arguments
- `filename`: file name. Assume it
    has one and only one kind of property group.
"""
macro getprop(filename)
    esc(quote
        h5open($(filename), "r") do f
            key = HDF5.names(f)[1] # only one property at top group
            d = read(f, key)
            args = [d[x] for x in Quaycle.prop_field_names[Symbol(key)]]
            eval(Expr(:call, Symbol(key), args...))
        end
    end)
end

"""
    store(filename::AbstractString, p::AbstractProperty)

Store property in HDF5.

## Arguments
- `filename::AbstractString`: file name to be used
- `p::AbstractProperty`: property to be saved
"""
function store(filename::AbstractString, p::AbstractProperty)
    h5open(filename, "w") do f
        g = g_create(f, description(p))
        for field in fieldnames(p)
            g[field] = getfield(p, Symbol(field))
        end
    end
end

"""
    @store filename::AbstractString p::AbstractProperty

Macro shortcut for storing property in HDF5.

## Arguments
- `filename::AbstractString`: file name to be used
- `p::AbstractProperty`: property to be saved
"""
macro store(filename, p)
    quote
        store($(esc(filename)), $(esc(p)))
    end
end
