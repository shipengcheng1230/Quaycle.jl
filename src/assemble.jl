## assemble the system derivative function

function assemble_ode(stype::Val{:okada}, fs::CentralSymmetryFS, sp::HomoFaultProperties, fp::FrictionalProperties; kwargs...)
    gf = greens_function(stype, fs.mesh, sp.λ, sp.μ, fs.dip, fs.faulttype; kwargs...)
    alloc = gen_alloc(stype, fs.mesh)
    dτ_dt! = gf_operator(stype, gf, alloc, sp.vpl)

    f! = (du, u, p, t) -> derivations!(du, u, p, tvar, se, fform)
end


@inline function dv_dθ_dt!(::CForm, dv::T, dθ::T, v::T, θ::T, se::StateEvolutionLaw) where {T<:AbstractVecOrMat}
    @fastmath @inbounds @simd for i = 1: prod(mp.dims)
        dμ_dθ = mp.σ[i] * mp.b[i] / θ[i]
        dμ_dv = mp.σ[i] * mp.a[i] / v[i]
        dθ[i] = dθ_dt(se, v[i], θ[i], mp.L[i])
        dv[i] = dv_dt(tvar.dτ_dt[i], dμ_dv, dμ_dθ, dθ[i], mp.η[i])
    end
end

function derivative_kernel!(du, u, p, t)

end
