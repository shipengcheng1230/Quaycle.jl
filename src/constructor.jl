function get_arg(name::Symbol, collection::Dict)
    x = get(collection, name, nothing)
    x == nothing && error("`$name` is not provided.")
    return x
end

broadcast_arg(x::Number, bcsize::NTuple) = x .* ones(bcsize...)
broadcast_arg(x::AbstractArray, bcsize::NTuple) = (size(x) == bcsize) ? x : error("$(size(x)) received, $bcsize required.")

parse_k(fa, gd, options::AbstractArray{<:Pair}) = parse_k(fa, gd, Dict(options))

function parse_k(fa, gd, options::Dict)
    if haskey(options, :filepath)
        filepath = options[:filepath]
        if isfile(filepath) && endswith(filepath, ".jld2")
            @info "Loading stiffness: $_k ..."
            @load filepath k
        else
            error("Invalid file path $filepath.")
        end
    elseif haskey(options, :μ)
        μ = options[:μ]
        λ = haskey(options, :λ) ? options[:λ] : μ
        ep = HomogeneousElasticProperties(λ=λ, μ=μ)
        @info "Calculating stiffness tensor ..."
        # `nrept` works only for 2D fault and is ignored by 1D case
        nrept = get(options, :nrept, 2)
        k = stiffness_tensor(fa, gd, ep; nrept=2)
    elseif haskey(options, :array)
        k = options[:array]
    else
        error("Insufficient options for stiffness tensor.")
    end
    return k
end

"""
    properties(fa::PlaneFaultDomain, gd::BoundaryElementGrid{dim}; _kwargs...) where {dim}

Establishing a material-properties-profile given by the fault domain and grids. User must provide the
    necessary parameters in according to the grid size specified or just a scalar for broadcasting.

## Arguments that are required:
- `a`: contrib from velocity.
- `b`: contrib from state.
- `L`: critical distance.
- `σ`: effective normal stress.
- `η`: radiation damping. It is recommended to set as ``μ / 2\\mathrm{Vs}`` where ``μ`` is *shear modulus* and ``\\mathrm{Vs}`` shear wave velocity.
- `vpl`: plate rate.
- `f0`: ref. frictional coeff.
- `v0`: ref. velocity.

## Arguments that need options
- `k`: stiffness tensor.

   (1) Providing *shear modulus* denoted as `μ` and *Lamé's first parameter* denoted as `λ` (same as `μ` if missing),
   then calculate it based on grid and fault domain, choosing parallel scheme if `nprocs() != 1`.
   (2) A valid file path to a *.jld2* that contains valid stiffness tensor. No verification will be performed here.
   (3) an `AbstractArray` represent the pre-calculated stiffness tensor. No verification will be performed here.
"""
properties(fa::PlaneFaultDomain{ftype, dim}, gd::BoundaryElementGrid, parameters::AbstractArray{<:Pair}) where {ftype<:PlaneFault, dim} = properties(fa, gd, Dict(parameters))
properties(;fault::PlaneFaultDomain{ftype, dim}, grid::BoundaryElementGrid, parameters::AbstractArray{<:Pair}) where {ftype<:PlaneFault, dim} = properties(fault, grid, parameters)

function properties(fa::PlaneFaultDomain{ftype, dim}, gd::BoundaryElementGrid{dim}, parameters::Dict) where {ftype<:PlaneFault, dim}
    bcsize = dim == 1 ? (gd.nξ,) : (gd.nx, gd.nξ)
    get_x = (name) -> broadcast_arg(get_arg(name, parameters), bcsize)
    a, b, L, σ, η = [get_x(name) for name in [:a, :b, :L, :σ, :η]]
    vpl, f0, v0 = [get_arg(name, parameters) for name in [:vpl, :f0, :v0]]
    k = parse_k(fa, gd, get_arg(:k, parameters))
    mp = PlaneMaterialProperties(dims=bcsize, a=a, b=b, L=L, k=k, σ=σ, η=η, vpl=vpl, f0=f0, v0=v0)
    @info "Fault material properties establised."
    return mp
end
