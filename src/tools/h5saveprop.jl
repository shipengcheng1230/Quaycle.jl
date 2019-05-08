export @read_prop, @save_prop, save_prop

macro read_prop(filename)
    esc(quote
        h5open($(filename), "r") do f
            c = read(f)
            ## only one property
            key = collect(keys(c))[1]
            d = c[key]
            args = [d[x] for x in fieldnames(Val(Symbol(key)))]
            eval(Expr(:call, Symbol(key), args...))
        end
    end)
end

function save_prop(filename::AbstractString, p::AbstractProperties)
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

macro save_prop(filename, p)
    esc(quote
        save_prop($(filename), $(p))
    end)
end
