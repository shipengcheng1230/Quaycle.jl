export @read_prop, @save_prop, save_prop

"Read property stored in HDF5."
macro read_prop(filename)
    esc(quote
        h5open($(filename), "r") do f
            key = names(f)[1] # only one property at top group
            d = read(f, key) 
            args = [d[x] for x in fieldnames(Val(Symbol(key)))]
            eval(Expr(:call, Symbol(key), args...))
        end
    end)
end

"Store property in HDF5."
function save_prop(filename::AbstractString, p::AbstractProperty)
    h5open(filename, "w") do f
        g = g_create(f, description(p))
        for field in fieldnames(p)
            g[field] = getfield(p, Symbol(field))
        end
    end
end

"Store property in HDF5."
macro save_prop(filename, p)
    quote
        save_prop($(esc(filename)), $(esc(p)))
    end
end
