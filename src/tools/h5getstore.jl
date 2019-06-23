export @getprop, @store, store

"Read property stored in HDF5."
macro getprop(filename)
    esc(quote
        h5open($(filename), "r") do f
            key = names(f)[1] # only one property at top group
            d = read(f, key)
            args = [d[x] for x in JuEQ.prop_field_names[Symbol(key)]]
            eval(Expr(:call, Symbol(key), args...))
        end
    end)
end

"Store property in HDF5."
function store(filename::AbstractString, p::AbstractProperty)
    h5open(filename, "w") do f
        g = g_create(f, description(p))
        for field in fieldnames(p)
            g[field] = getfield(p, Symbol(field))
        end
    end
end

"Store property in HDF5."
macro store(filename, p)
    quote
        store($(esc(filename)), $(esc(p)))
    end
end
