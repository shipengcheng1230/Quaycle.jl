## plastic deformation
export DislocationCreep, DiffusionCreep, Peierls

abstract type DeformationMechanism end

# Karato 2008, Deformation of Earth Materials
struct DiffusionCreep <: DeformationMechanism end
struct DislocationCreep <: DeformationMechanism end
struct Peierls <: DeformationMechanism end # high stress / low temperature

const 𝙍 = float(PhysicalConstants.CODATA2014.R).val # ideal gas constant

# stress driven plastic deformation
function dϵ_dt(::DislocationCreep, A, σ, τ, n, fH₂0, r, α, ϕ, Q, P, Ω, T)
    A * τ^(n-1) * σ * fH₂0^r * exp(α * ϕ) * exp(-(Q + P * Ω) / 𝙍 / T)
end

function dϵ_dt(::DiffusionCreep, A, σ, d, m, fH₂0, r, α, ϕ, Q, P, Ω, T)
    A * σ * d^(-m) * fH₂0^r * exp(α * ϕ) * exp(-(Q + P * Ω) / 𝙍 / T)
end
