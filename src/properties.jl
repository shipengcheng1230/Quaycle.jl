## properties interface

export init_fri_prop, init_sys_prop, read_properties, save_properties

import Base.fieldnames

abstract type AbstractProperties end

@with_kw struct FrictionalProperties{U<:AbstractVecOrMat} <: AbstractProperties
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    σ::U # effective normal stress
end


@with_kw struct SystemProperties{T<:Number} <: AbstractProperties
    η::T # radiation damping
    vpl::T # plate rate, unlike pure Rate-State Friction simulation, here is restrained to be constant
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    @assert η > 0
    @assert vpl > 0
    @assert f0 > 0
    @assert v0 > 0
end


fieldnames(p::FrictionalProperties) = ("a", "b", "L", "σ")

fieldnames(p::SystemProperties) = ("η", "vpl", "f0", "v0")

description(p::FrictionalProperties) = "friction"

description(p::SystemProperties) = "system"

init_sys_prop(η, vpl, f0=0.6, v0=1e-6) = SystemProperties(promote(η, vpl, f0, v0)...)

init_fri_prop(mesh::SimpleLineGrid) = FrictionalProperties([Vector{eltype(mesh.Δξ)}(undef, mesh.nξ) for _ in 1: 4]...)

init_fri_prop(mesh::SimpleRectGrid) = FrictionalProperties([Matrix{eltype(mesh.Δx)}(undef, mesh.nx, mesh.nξ) for _ in 1: 4]...)


function read_properties(filepath::AbstractString)
    c = h5open(filepath, "r") do f
        read(f)
    end
    props = Dict()

    if haskey(c, "system")
        d = c["system"]
        props["system"] = SystemProperties(d["η"], d["vpl"], d["f0"], d["v0"])
    end

    if haskey(c, "friction")
        d = c["friction"]
        props["friction"] = FrictionalProperties(d["a"], d["b"], d["L"], d["σ"])
    end

    return props
end


function save_properties(filename::AbstractString, p::AbstractProperties)
    h5open(filename, "cw") do f
        c = read(f)
        if haskey(c, description(p))
            g = c[description(p)]
        else
            g = g_create(f, description(p))
        end
        for f in fieldnames(p)
            g[f] = getfield(p, Symbol(f))
        end
    end
end


function save_properties(filename::AbstractString, ps::Vector{AbstractProperties})
    for p in ps
        save_properties(filename, p)
    end
end
