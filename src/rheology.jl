abstract type DeformationMechanism end

# Karato 2008, Deformation of Earth Materials
struct DiffusionCreep <: DeformationMechanism end
struct DislocationCreep <: DeformationMechanism end
struct Peierls <: DeformationMechanism end # high stress / low temperature

const _R = float(CODATA2014.R).val # ideal gas constant

# stress driven plastic deformation
function dϵ_dt(::DiffusionCreep, A, σ, τ, n, fH₂0, r, Q, T)
    A * τ^(n-1) * σ * fH₂0^r * exp(-Q / _R / T)
end

function dϵ_dt(::DislocationCreep, A, σ, d, m, fH₂0, r, Q, T)
    A * σ * d^(-m) * fH₂0^r * exp(-Q / _R / T)
end
