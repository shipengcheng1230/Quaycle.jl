# This is a demo for SEAS BP1, inspired by Yajing Liu's code
# 2D vertical anti-shear-plain

## let's do parallel
addprocs(4)
# use dc3d
dc3d_file = joinpath(dirname(@__DIR__), "src/dc3d.jl")
# need string insertion here to work around
@everywhere include($dc3d_file)
using DifferentialEquations
using Plots

## parameters are set using BP1 requirements
const ρ = 2670.0
const g = 9.8
const cs = 11.08e7
const σ = 500.0
const a0 = 0.010
const amax = 0.025
const b0 = 0.015
const Dc = 8.0
const Vp = 31.53
const Vinit = 31.53
const V0 = 31.53e3
const f0 = 0.6
const H = 15.0
const h = 3.0
const Wf = 40.0
const Δz = 25.0e-3
const tf = 3000.0

# parameters implicit by above
const μ = 0.32
const λ = μ
const α = (λ + μ) / (λ + 2μ)
const η = μ / 2cs

# some settings
const ngrid = Int(Wf / Δz)
const depth = 0.
const al = [-500., 500.]
const disl = [1., 0., 0.]
const dip = 90.0
const cdip = cospi(dip / 180)
const sdip = sinpi(dip / 180)

## stiff matrix, parallelly assembled
K = SharedArray{Float64}(ngrid, ngrid)

# row: obs depths; col: fault patches
function fill_stiff_matrix!(K)
    for i = 1: ngrid
        aw = [-Δz * i, -Δz * (i - 1)]
        # if no reduction op, @parallel will not wait for finish without @sync
        @sync @parallel for j = 1: ngrid
            x = 0.
            y = -(j - 0.5) * Δz * cdip
            z = -(j - 0.5) * Δz * sdip
            u = dc3d_okada(x, y, z, α, depth, dip, al, aw, disl)
            tn = μ * (u[5] + u[7])
            K[j, i] = tn / disl[1]
        end
    end
end

@time fill_stiff_matrix!(K)
## build depth profile
z = (collect(1: ngrid) - 0.5) * Δz
az = fill!(zeros(z), a0)
az[z .≥ (H + h)] = amax
az[H .< z .< H + h] = a0 + (amax - a0) / (h / Δz)* collect(1: Int(h / Δz))

# initial condition
δz = zeros(z)
τ0 = σ * amax * asinh(Vinit / 2V0 * exp((f0 + b0 * log(V0 / Vinit)) / amax)) + η * Vinit
τz = fill!(zeros(z), τ0)
θz = @. Dc / V0 * exp(az / b0 * log(2V0 / Vinit * sinh((τz - η * Vinit) / az / σ)) - f0 / b0)
vz = fill!(zeros(z), Vinit)

## derivation
const K2 = convert(Array, K)

@views function rsf_odes!(du, u, p, t)
    kk, vp, aa, bb, dc, eta = p
    dμ_dt = kk * (vp .- u[:, 1])
    dθ_dt = @. 1 - u[:, 1] * u[:, 2] / dc
    dμ_dθ = @. bb / u[:, 2]
    dμ_dv = @. aa / u[:, 1]
    dv_dt = @. (dμ_dt - dμ_dθ * dθ_dt) / (dμ_dv + eta)
    du[:, 1] = dv_dt
    du[:, 2] = dθ_dt
end

function rsf_odes_regularized!(du, u, p, t)
    ϕ1 = @. log(V0 * u[:, 2] / Dc)
    ϕ2 = @. u[:, 1] / 2V0
    ϕ3 = @. (f0 + b0 * ϕ1) / az
    ϕ4 = @. sqrt(1 + (ϕ2 * exp(ϕ3))^2)
    dμ_dv = @. az / ϕ4 * exp(ϕ3) / 2V0
    dμ_dθ = @. b0 / u[:, 2] * ϕ2 * exp(ϕ3) / ϕ4
    dθ_dt = @. 1 - u[:, 1] * u[:, 2] / Dc
    dμ_dt = K * (Vp - u[:, 1])
    dv_dt = @. (dμ_dt - dμ_dθ * dθ_dt) / (dμ_dv + η)
    du[:, 1] = dv_dt
    du[:, 2] = dθ_dt
end

u0 = hcat(vz, θz)
tspan = (0., 1.0)
prob = ODEProblem(rsf_odes!, u0, tspan, (K2, Vp, az, b0, Dc, η))
@time sol = solve(prob, Tsit5(), reltol=1e-3, abstol=1e-3)

##
plot(sol.t, sol[1, 1, :])
plot(az, z)
