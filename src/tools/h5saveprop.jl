export @read_prop, @save_prop, save_prop

"Read property stored in HDF5."
macro read_prop(filename)
    quote
        h5open($(esc(filename)), "r") do f
            c = read(f)
            ## only one property
            key = collect(keys(c))[1]
            d = c[key]
            args = [d[x] for x in fieldnames(Val(Symbol(key)))]
            eval(Expr(:call, Symbol(key), args...))
        end
    end
end

"Store property in HDF5."
function save_prop(filename::AbstractString, p::AbstractProperty)
    h5open(filename, "w") do f
        c = read(f)
        if haskey(c, description(p))
            g = c[description(p)]
        else
            g = g_create(f, description(p))
        end
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
