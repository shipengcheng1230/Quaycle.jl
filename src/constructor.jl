"""
    properties(fa::PlaneFaultDomain, gd::BoundaryElementGrid{dim}; _kwargs...) where {dim}

Establishing a material-properties-profile given by the fault domain and grids. User must provide the
    necessary parameters in according to the grid size specified or just a scalar for broadcasting.

## Arguments that are needed:
- `a`: contrib from velocity.
- `b`: contrib from state.
- `L`: critical distance.
- `σ`: effective normal stress.
- `vpl`: plate rate.
- `f0`: ref. frictional coeff.
- `v0`: ref. velocity.

## Arguments that are optional
- `k`: stiffness tensor. If `:auto`, it will automatically calculate by seeking `λ` and `μ` otherwise should be a valid file path to a `*.jld2` or an `AbstractArray`.
- `η`: radiation damping. If `:auto`, it will automatically seek `μ` and `vs` and use ``μ / 2\\mathrm{Vs}``.
- `vs`: shear wave velocity.
- `λ`: Lamé's first parameter
- `μ`: shear modulus
"""
function properties(fa::PlaneFaultDomain, gd::BoundaryElementGrid{dim}; _kwargs...) where {dim}

    function get_k()
        _k = get(kwargs, :k, :auto)
        _tk = typeof(_k)
        if _tk <: AbstractString && isfile(_k) && endswith(_k, ".jld2")
            @info "Loading stiffness: $_k ..."
            @load _k k
        elseif _tk <: AbstractArray
            k = _k
        elseif _k == :auto
            λ = args_get_expand(:λ, kwargs, (), false)
            μ = args_get_expand(:μ, kwargs, (), false)
            ep = HomogeneousElasticProperties(λ=μ, μ=μ)
            @info "Calculating stiffness tensor..."
            k = stiffness_tensor(fa, gd, ep)
        else
            error("""
                Invalid option: $_k, should be:
                (1) A valid file path to a `*.jld2` file
                (2) An `AbstractArray`
                (3) Symbol `:auto`.
                (4) NamedTuple containts `:nrept`.
            """)
        end

        path = get(kwargs, :savek, nothing)
        if typeof(path) <: AbstractString && ispath(path)
            @info "Saving stiffness: $path ..."
            @save path k
        end
        return k
    end

    function get_η()
        _η = get(kwargs, :η, nothing)
        if _η == :auto || _η == nothing
            μ = args_get_expand(:μ, kwargs, (), false)
            vs = args_get_expand(:vs, kwargs, gsize)
            η = μ ./ 2vs
        else
            η = args_get_expand(:η, kwargs, gsize)
        end
    end

    gsize = dim == 1 ? (gd.nξ,) : (gd.nx, gd.nξ)
    xind = dim == 1 ? fill(true, gd.nξ) : gd.xind
    kwargs = _kwargs.data
    a = args_get_expand(:a, kwargs, gsize)
    b = args_get_expand(:b, kwargs, gsize)
    L = args_get_expand(:L, kwargs, gsize)
    σ = args_get_expand(:σ, kwargs, gsize)

    vpl = args_get_expand(:vpl, kwargs, (), false)
    f0 = args_get_expand(:f0, kwargs, (), false)
    v0 = args_get_expand(:v0, kwargs, (), false)

    k = get_k()
    η = get_η()

    @info "Establishing material properties..."
    mp = PlaneMaterialProperties(comsize=gsize ,a=a, b=b, L=L, k=k, σ=σ, η=η, vpl=vpl, f0=f0, v0=v0)
    return mp
end
