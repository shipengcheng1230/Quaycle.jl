## properties interface

export init_friction_prop, init_fault_prop, read_properties, save_properties

import Base.fieldnames

abstract type AbstractProperties end

@with_kw struct FrictionalProperties{
    U<:AbstractVecOrMat, FForm<:FrictionLawForm, SEL<:StateEvolutionLaw
    } <: AbstractProperties
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    σ::U # effective normal stress
    flf::FForm = RForm() # frictional law format
    sel::SEL = DieterichStateLaw() # state evolution law
end


@with_kw struct HomoFaultProperties{T<:Number} <: AbstractProperties
    λ::T # lame first para
    μ::T # shear modulus
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

ff = RForm()
ex = Expr(:call, :RForm)
ff2 = eval(Expr(:call, Symbol("DieterichStateLaw")))
eval((:call, :RForm))


fieldnames(p::FrictionalProperties) = ("a", "b", "L", "σ", "flf", "sel")
fieldnames(p::HomoFaultProperties) = ("λ", "μ", "η", "vpl", "f0", "v0")

description(p::FrictionalProperties) = "friction"
description(p::HomoFaultProperties) = "faultspace"

init_fault_prop(η, vpl, f0=0.6, v0=1e-6) = HomoFaultProperties(promote(λ, μ, η, vpl, f0, v0)...)
init_friction_prop(mesh::SimpleLineGrid) = FrictionalProperties([Vector{eltype(mesh.Δξ)}(undef, mesh.nξ) for _ in 1: 4]...)
init_friction_prop(mesh::SimpleRectGrid) = FrictionalProperties([Matrix{eltype(mesh.Δx)}(undef, mesh.nx, mesh.nξ) for _ in 1: 4]...)

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
        props["friction"] = FrictionalProperties(
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
