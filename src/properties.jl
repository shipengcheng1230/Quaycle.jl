## properties interface

export init_friction_prop, init_fault_prop, read_properties, save_properties
export HomoFaultProperties, RSFrictionalProperties
export description

import Base.fieldnames

abstract type AbstractProperties end

@with_kw struct RSFrictionalProperties{
    U<:AbstractVecOrMat, FForm<:FrictionLawForm, SEL<:StateEvolutionLaw
    } <: AbstractProperties
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    σ::U # effective normal stress
    flf::FForm = RForm() # frictional law format
    sel::SEL = DieterichStateLaw() # state evolution law

    @assert size(a) == size(b)
    @assert size(b) == size(L)
    @assert size(L) == size(σ)
end


@with_kw struct HomoFaultProperties{T<:Number} <: AbstractProperties
    λ::T # Lamé first constants
    μ::T # Lamé second constants
    η::T # radiation damping
    vpl::T # plate rate, unlike pure Rate-State Friction simulation, here is restrained to be constant
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    @assert λ > 0
    @assert μ > 0
    @assert η > 0
    @assert vpl > 0
    @assert f0 > 0
    @assert v0 > 0
end

fieldnames(p::RSFrictionalProperties) = ("a", "b", "L", "σ", "flf", "sel")
fieldnames(p::HomoFaultProperties) = ("λ", "μ", "η", "vpl", "f0", "v0")

description(p::RSFrictionalProperties) = "friction"
description(p::HomoFaultProperties) = "faultspace"

init_fault_prop(λ, μ, η, vpl, f0=0.6, v0=1e-6) = HomoFaultProperties(promote(λ, μ, η, vpl, f0, v0)...)

init_friction_prop(mesh::SimpleLineGrid) = RSFrictionalProperties(
    [Vector{eltype(mesh.Δξ)}(undef, mesh.nξ) for _ in 1: 4]...,
    RForm(), DieterichStateLaw())
init_friction_prop(mesh::SimpleRectGrid) = RSFrictionalProperties(
    [Matrix{eltype(mesh.Δx)}(undef, mesh.nx, mesh.nξ) for _ in 1: 4]...,
    RForm(), DieterichStateLaw())
init_friction_prop(fs::CentralSymmetryFS) = init_friction_prop(fs.mesh)

function read_properties(filepath::AbstractString)
    c = h5open(filepath, "r") do f
        read(f)
    end
    props = Dict()

    if haskey(c, "faultspace")
        d = c["faultspace"]
        props["faultspace"] = HomoFaultProperties(d["λ"], d["μ"], d["η"], d["vpl"], d["f0"], d["v0"])
    end

    if haskey(c, "friction")
        d = c["friction"]
        props["friction"] = RSFrictionalProperties(
            d["a"], d["b"], d["L"], d["σ"],
            eval(Expr(:call, Symbol(d["flf"]))),
            eval(Expr(:call, Symbol(d["sel"]))),
            )
    end

    return props
end

function save_properties(filename::AbstractString, p::AbstractProperties; option="cw")
    h5open(filename, option) do f
        c = read(f)
        if haskey(c, description(p))
            g = c[description(p)]
        else
            g = g_create(f, description(p))
        end
        for field in fieldnames(p)
            g[field] = getfield(p, Symbol(field)) |> __h5_compatible_converter__
        end
    end
end

__h5_compatible_converter__(x) = x
# h5 cannot write substring
__h5_compatible_converter__(x::FrictionLawForm) = split(string(typeof(x)), '.')[end] |> string
__h5_compatible_converter__(x::StateEvolutionLaw) = split(string(typeof(x)), '.')[end] |> string

function save_properties(filename::AbstractString, ps::Vector{AbstractProperties})
    for p in ps
        save_properties(filename, p)
    end
end
