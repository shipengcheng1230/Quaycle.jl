# !!! note
#     This example is from [Benchmark Problem 1](https://www.scec.org/publication/8214) (hence referred as BP1).
#
# ### Define Parameters

# First, we load the package
using JuEQ
using Plots

# Instead of using SI unit, we refactor ours into the follow:

ms2mmyr = 365 * 86400 * 1e3 # convert velocity from m/s to mm/yr
ρ = 2670.0 # density [kg/m³]
vs = 3464.0 # shear wave velocity [m/s]
σ = 500.0 # effective normal stress [bar]
a0 = 0.010 # frictional paramter `a` in vw zone
amax = 0.025 # frictional paramter `a` in vs zone
b0 = 0.015 # frictional paramter `b`
L = 8.0 # critical distance [mm]
vpl = 1e-9 * ms2mmyr # plate rate [mm/yr]
vinit = 1e-9 * ms2mmyr # initial velocity [mm/yr]
v0 = 1e-6 * ms2mmyr # reference velocity [mm/yr]
f0 = 0.6 # reference frictional coefficient
H = 15.0 # vw zone [km]
h = 3.0 # vw-vs changing zone [km]
Wf = 40.0 # fault depth [km]
Δz = 100.0e-3 # grid size interval [km]
tf = 400.0; # simulation time [yr]

# !!! warning
#     Make sure your units are consistent across the whole variable space. Pontenial imporvement may
#     incoporate [Unitful.jl](https://github.com/ajkeller34/Unitful.jl) package.

# Then we arrive at some parameters that are implicit by above:

μ = vs^2 * ρ / 1e5 / 1e6 # shear modulus [bar·km/mm]
λ = μ # poisson material
η = μ / 2(vs * 1e-3 * 365 * 86400)
ngrid = round(Int, Wf / Δz); # number of grids

# Now, we start to construct our model using parameters above. First, we create a 'fault' by specifying fault type and depth:
# !!! tip
#     Here, we do not need to provide `dip` for strike-slip fault as it automatically choose `90`. See [`fault`](@ref).

# ### Construct Model

fa = fault(StrikeSlipFault, Wf);

# Next, we generate the grid regarding the `fault` we just created by giving number of grids:
# !!! note
#     This package use `ξ` for denoting downdip coordinate and `x` for along-strike one. See [`discretize`](@ref).

gd = discretize(fa; nξ=ngrid);

# Next, we construct the required frictional parameter profile:

z = -gd.ξ
az = fill(a0, size(z))
az[z .≥ (H + h)] .= amax
az[H .< z .< H + h] = a0 .+ (amax - a0) / (h / Δz) * collect(1: Int(h / Δz));

# Then, we provide the required initial condition satisfying uniform slip distribution over the depth:

τ0 = σ * amax * asinh(vinit / 2v0 * exp((f0 + b0 * log(v0 / vinit)) / amax)) + η * vinit
τz = fill(τ0, size(z))
θz = @. L / v0 * exp(az / b0 * log(2v0 / vinit * sinh((τz - η * vinit) / az / σ)) - f0 / b0)
vz = fill(vinit, size(z))
u0 = hcat(vz, θz);

# Let's simulate only the first 200 years:

tspan = (0., 200.);

# Finally, we provide the material properties w.r.t. our 'fault', 'grid' as well as other necessary parameters predefined
# using the same grid size & dimension:

mp = properties(;fault=fa, grid=gd, parameters=[:a=>az, :b=>b0, :L=>L, :σ=>σ, :η=>η, :k=>[:λ=>λ, :μ=>μ], :vpl=>vpl, :f0=>f0, :v0=>v0]);

# !!! tip
#     Check [`properties`](@ref) for extended options.

# Check our profile now:

plot([mp.a, mp.b], z, label=["a", "b"], yflip=true, ylabel="Depth (km)")

# We then contruct the `ODEProblem` as following by stating which state evolution law to use and frcitonal law form,
# plus initial condition and simulation time:

prob = EarthquakeCycleProblem(gd, mp, u0, tspan; se=DieterichStateLaw(), fform=RForm());

# ### Solve Model
# We then solve the ODEs:

sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6);

# !!! tip
#     For details of solving options, see [here](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html).

# !!! tip
#     Raise the accuracy option if you get instability when solving these ODEs.

# ### Results
# The first event happens at around 196 year:

maxv = max_velocity(sol)
plot(sol.t, log10.(maxv / ms2mmyr), xlabel="Time (year)", ylabel="Max Velocity (log10 (m/s))", xlims=(190, 200), label="")

# !!! note
#     Click [here](https://plot.ly/~shipengcheng/3/depth-km-vs-slip-m/) for the slip evolution over 3000 years simulation.
#     It may need some time to load the page.
