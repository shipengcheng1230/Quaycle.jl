# !!! note
#     This example is from [Benchmark Problem 1](https://www.scec.org/publication/8214) (hence referred as BP1).
#

# First, we load the package
using JuEQ
using Plots

# The prerequisite parameters in this benchmark are list below:

ms2mmyr = 365 * 86400 * 1e3 # convert velocity from m/s to mm/yr
ρ = 2670.0 # density [kg/m³]
vs = 3464.0 # shear wave velocity [m/s]
σ0 = 500.0 # effective normal stress [bar]
a0 = 0.010 # frictional paramter `a` in vw zone
amax = 0.025 # frictional paramter `a` in vs zone
b0 = 0.015 # frictional paramter `b`
L0 = 8.0 # critical distance [mm]
vpl = 1e-9 * ms2mmyr # plate rate [mm/yr]
vinit = 1e-9 * ms2mmyr # initial velocity [mm/yr]
v0 = 1e-6 * ms2mmyr # reference velocity [mm/yr]
f0 = 0.6 # reference frictional coefficient
H = 15.0 # vw zone [km]
h = 3.0 # vw-vs changing zone [km]
Wf = 40.0 # fault depth [km]
Δz = 100.0e-3 # grid size interval [km]
tf = 200.0; # simulation time [yr]
nothing

# !!! warning
#     Make sure your units are consistent across the whole variable space.

# Then we arrive at some parameters that are implicit by above:

μ = vs^2 * ρ / 1e5 / 1e6 # shear modulus [bar·km/mm]
λ = μ # poisson material
η = μ / 2(vs * 1e-3 * 365 * 86400)
ngrid = round(Int, Wf / Δz); # number of grids
nothing

# First, create a **fault space** by specifying fault type, depth
# and with the desired discretization interval.
# !!! tip
#     Here, we do not need to provide `dip` for strike-slip fault as it automatically choose `90`. See [`fault`](@ref).

fa = fault(Val(:LineOkada), STRIKING(), 40.0, Δz)
nothing

# Then, provide the material properties w.r.t. our 'fault space'.
a = a0 .* ones(fa.mesh.nξ)
a[-fa.mesh.z .≥ (H + h)] .= amax
a[H .< -fa.mesh.z .< H + h] .= a0 .+ (amax - a0) / (h / Δz) * collect(1: Int(h / Δz))
b = b0 .* ones(fa.mesh.nξ)
L = L0 .* ones(fa.mesh.nξ)
σ = σ0 .* ones(fa.mesh.nξ)
prop = ElasticRSFProperty(a=a, b=b, L=L, σ=σ, λ=λ, μ=μ, vpl=vpl, f0=f0, v0=v0, η=η)
nothing

# Next, construct the initial condition and ODE problem using Okada's Green's function.
τ0 = σ0 * amax * asinh(vinit / 2v0 * exp((f0 + b0 * log(v0 / vinit)) / amax)) + η * vinit
τz = fill(τ0, size(fa.mesh.z))
θz = @. L / v0 * exp(a / b0 * log(2v0 / vinit * sinh((τz - η * vinit) / a / σ)) - f0 / b0)
vz = fill(vinit, size(fa.mesh.ξ))
u0 = ArrayPartition(vz, θz)
prob, = assemble(fa, prop,  u0, (0.0, tf))
nothing

# Check our depth profile now.

plot([a, b], fa.mesh.z, label=["a", "b"], yflip=true, ylabel="Depth (km)")

# Afterwards, solve ODE thanks to [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)

sol = solve(prob, TsitPap8(), reltol=1e-6, abstol=1e-6)
nothing

# !!! tip
#     Raise the accuracy option or switch to other algorithms if you get instability when solving these ODEs.

# Finally, check the results. The first event happens at around 196 year:

maxv = max_velocity(sol)
plot(sol.t, log10.(maxv / ms2mmyr), xlabel="Time (year)", ylabel="Max Velocity (log10 (m/s))", xlims=(190, 200), label="")

# !!! note
#     Click [here](https://plot.ly/~shipengcheng/3/depth-km-vs-slip-m/) for the slip evolution over 3000 years simulation.
#     It may need some time to load the page.
