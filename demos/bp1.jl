# This is a demo for SEAS BP1, inspired by Yajing Liu's code
# 2D vertical anti-shear-plain

## let's do parallel
ncores = 4
addprocs(ncores)

## parameters are set using BP1 requirements
const ms2mmyr = 365 * 86400 * 1e3
const ρ = 2670.0 # kg/m³
const cs = 3464.0 # m/s
const σ = 500.0 # Bar
const a0 = 0.010
const amax = 0.025
const b0 = 0.015
const Dc = 8.0 # mm
const Vp = 1e-9 * ms2mmyr # mm/yr
const Vinit = 1e-9 * ms2mmyr # mm/yr
const V0 = 1e-6 * ms2mmyr # mm/yr
const f0 = 0.6
const H = 15.0 # km
const h = 3.0 # km
const Wf = 40.0 # km
const Δz = 25.0e-3 # km
const tf = 3000.0 # yr

# parameters implicit by above
const μ = cs^2 * ρ / 1e5 / 1e6 # Bar·km/mm
const λ = μ # poisson material
const α = (λ + μ) / (λ + 2μ)
const η = μ / 2(cs * 1e-3 * 365 * 86400) # Bar·yr/mm
const ngrid = Int(Wf / Δz)

## calculate stiffness matrix
dc3d_file = joinpath(dirname(@__DIR__), "src/dc3d.jl")
# need string insertion here to work around
@everywhere include($dc3d_file)

K = SharedArray{Float64}(ngrid, ngrid)

# fault settings
const depth = 0.0
const al = [-2000.0, 2000.] # periodic boundary condition
const disl = [1., 0., 0.]
const dip = 90.0
const cdip = cospi(dip / 180)
const sdip = sinpi(dip / 180)

# row: obs along depth; col: fault patches
function fill_stiffness_matrix!(K)
    for i = 1: ngrid
        aw = [-Δz * i, -Δz * (i - 1)]
        # if no reduction op, @parallel will not wait for finish without @sync
        @sync @parallel for j = 1: ngrid
            x = 0.
            y = -(j - 0.5) * Δz * cdip
            z = -(j - 0.5) * Δz * sdip
            u = dc3d_okada(x, y, z, α, depth, dip, al, aw, disl)
            tn = μ * (u[5] + u[7])
            # mode 3 crack
            K[j, i] = tn / disl[1]
        end
    end
end

@time fill_stiffness_matrix!(K)
K = sdata(K)

# don't have parallel ODE solver for now
rmprocs(ncores)

## save or load your stiffness matrix
using JLD2, FileIO
@save joinpath(@__DIR__, "bp1_stiff.jld2") K
# load existing solution
@load joinpath(@__DIR__, "bp1_stiff.jld2") K

## profile settings
z = (collect(1: ngrid) - 0.5) * Δz
az = fill!(zeros(z), a0)
az[z .≥ (H + h)] = amax
az[H .< z .< H + h] = a0 + (amax - a0) / (h / Δz) * collect(1: Int(h / Δz))
# global constants are more efficient in case of frequent access
const a = az

# initial condition
δz = zeros(z)
τ0 = σ * amax * asinh(Vinit / 2V0 * exp((f0 + b0 * log(V0 / Vinit)) / amax)) + η * Vinit
τz = fill!(zeros(z), τ0)
θz = @. Dc / V0 * exp(az / b0 * log(2V0 / Vinit * sinh((τz - η * Vinit) / az / σ)) - f0 / b0)
vz = fill!(zeros(z), Vinit)

## build ODEs
using DifferentialEquations

ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv = [zeros(Float64, ngrid) for _ in 1: 5]

function f_full!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv)
    v = @view u[:, 1]
    θ = @view u[:, 2]
    dv = @view du[:, 1]
    dθ = @view du[:, 2]
    dδ = @view du[:, 3]

    # make sure θ don't go below 0.0
    clamp!(θ, 0.0, Inf)

    # using regularized form
    @. ϕ1 = exp((f0 + b0 * log(V0 * θ / Dc)) / a) / 2V0
    @. ϕ2 = σ * ϕ1 / sqrt(1 + (v * ϕ1)^2)
    @. dμ_dv = a * ϕ2
    @. dμ_dθ = b0 / θ * v * ϕ2
    @. dθ = 1 - v * θ / Dc
    A_mul_B!(dμ_dt, K, Vp - v)
    @. dv = (dμ_dt - dμ_dθ * dθ) / (dμ_dv + η)
    @. dδ = v
end

# using closures to minimize allocations
f! = (du, u, p, t) -> f_full!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv)

u0 = hcat(vz, θz, δz)
tspan = (0., tf)
prob = ODEProblem(f!, u0, tspan)

# solution saving options, note that all solutions are stored at memory at runtime for now
const Δt_seismic = 0.1 / 365 / 86400 # 0.1 s
const Δt_aseismic = 0.1 # 0.1 year
const seismic_slip_rate = 1e-3 * ms2mmyr # above which an event is deemed as seismic

# Note that by this `affect!` function, seismic solutions are not rigorously intepolated at predefined Δt
function affect!(integrator::DEIntegrator)
    t_saved = integrator.sol.t[end]
    t = integrator.t
    if t - t_saved ≥ Δt_seismic
        add_saveat!(integrator, t)
    end
end

# identify seismic period
function condition(u, t, integrator)
    maxv = maximum(u[:, 1])
    maxv >= seismic_slip_rate
end

cb = DiscreteCallback(condition, affect!, save_positions=(false, false))

# this is the only kind of algorithms that works on this problem efficiently
@time sol = solve(prob, OwrenZen5(), reltol=1e-6, abstol=1e-6, saveat=Δt_aseismic, callback=cb)

## save our solution for plotting
using HDF5
using RecursiveArrayTools

# convert higher dimension array to matrix
uu = VectorOfArray(sol.u)
u = convert(Array, uu)

# regularized form
function shear_stress(v, θ, a)
    @. σ * a * asinh(v / 2V0 * exp((f0 + b0 * log(V0 * θ / Dc)) / a))
end

# organize various outputs
v = @view u[:, 1, :]
θ = @view u[:, 2, :]
s = @view u[:, 3, :]
τ = shear_stress(v, θ, a)

v /= ms2mmyr
θ *= 365 * 86400
s /= 1e3
τ /= 10.

# write results
h5open(joinpath(@__DIR__, "bp1_solution.h5"), "w") do f
    g = g_create(f, "bp1")
    g["t"] = sol.t # yr
    g["velocity"] = v # m/s
    g["state"] = θ # s
    g["shear_stress"] = τ # MPa
    g["slip"] = s # m
end
