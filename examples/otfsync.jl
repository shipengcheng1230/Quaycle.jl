# !!! note
#     This example is adapted from [Wei, 2016 AGU](https://agu.confex.com/agu/fm16/meetingapp.cgi/Paper/132124)

# !!! tip
#     It will automatically use parallel scheme if `nprocs() ≂̸ 1` when building stiffness tensor. To do so:
#     ```julia
#     using Distributed
#     addprocs(4); # add # of cores you desire
#     using JuEQ
#     ```

# First, we load the package and define some basic parameters:

using JuEQ
using Plots

ms2mmyr = 365 * 86400 * 1e3
ρ = 2670.0 # kg/m³
cs = 3044.0 # m/s
vpl = 100.0 # mm/yr
v0 = 3.2e4 # mm/yr
f0 = 0.6;

# Then we come to parameters implicit by above:

μ = 0.3 # Bar·km/mm
λ = μ # poisson material
α = (λ + μ) / (λ + 2μ)
η = μ / 2(cs * 1e-3 * 365 * 86400); # Bar·yr/mm

# Create a fault:

fa = fault(StrikeSlipFault, (80., 10.));

# Generate grids:

gd = discretize(fa; nx=160, nξ=20, buffer=0.);

# !!! tip
#     It is recommended (from Yajing Liu's personal communication) to add buffer zones adjacent the horizontal edges
#     to immitate *zero* dislocation at the ridge region.
#     Basically, it affects how the stiffness tensor are periodically summed. To what extent it alters the results remains further testing.

#     Under such circumstance, the grid will expand outside the actual fault domain. However, in the buffer zone, relative
#     velocity is *zero* at all time, so are the derivatives of velocity & state variable. Users may need to extract the solutions
#     in the actually fault domain using `xs_index::BitArray`, one of the fields of the grid structure
#     that indicate which part of along-strike is the actually fault domain.

# Time for us to establish frictional parameters profile:

a = 0.015 .* ones(gd.nx, gd.nξ)
b = 0.0115 .* ones(gd.nx, gd.nξ)
left_patch = @. -25. ≤ gd.x ≤ -5.
right_patch = @. 5. ≤ gd.x ≤ 25.
vert_patch = @. -6. ≤ gd.z ≤ -1.
b[xor.(left_patch, right_patch), vert_patch] .= 0.0185
amb = a - b
σmax = 500.
σ = [min(σmax, 15. + 180. * z) for z in -gd.z]
σ = Matrix(repeat(σ, 1, gd.nx)')
L = 12.;

# Check our profile:

p1 = heatmap(amb',
    xticks=(0: 10/gd.Δx: gd.nx, -fa.span[1]/2: 10: fa.span[1]/2),
    yticks=(0: 5/gd.Δξ: gd.nξ, 0: -5: -fa.span[2]),
    yflip=true, color=:isolum, aspect_ratio=2, title="a-b"
    );

p2 = heatmap(σ',
    xticks=(0: 10/gd.Δx: gd.nx, -fa.span[1]/2: 10: fa.span[1]/2),
    yticks=(0: 5/gd.Δξ: gd.nξ, 0: -5: -fa.span[2]),
    yflip=true, color=:isolum, aspect_ratio=2, title="\\sigma"
    );

plot(p1, p2, layout=(2, 1))

# Construct our material property profile:

mp = properties(fa, gd; a=a, b=b, L=L, k=:auto, σ=σ, η=η, v0=v0, f0=f0, vpl=vpl, λ=λ, μ=μ);

# Provide the initial condition:

vinit = vpl .* ones(gd.nx, gd.nξ)
θ0 = L ./ vinit ./ 1.1
u0 = cat(vinit, θ0, dims=3);

# Get our ODEs problem:

prob = EarthquakeCycleProblem(gd, mp, u0, (0., 18.); se=DieterichStateLaw(), fform=CForm());

# Solve the model:

@time sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6);

# Take a look at the max velocity:

maxv = max_velocity(sol)
plot(sol.t, log10.(maxv / ms2mmyr), xlabel="Time (year)", ylabel="Max Velocity (log10 (m/s))", label="")

# View some snapshots to see the rupture (quasi-dynamic) patterns:

ind = argmax(maxv)
myplot = (ind) -> heatmap(log10.(sol.u[ind][:,:,1]./ms2mmyr)',
    xticks=(0: 10/gd.Δx: gd.nx, -fa.span[1]/2: 10: fa.span[1]/2),
    yticks=(0: 5/gd.Δξ: gd.nξ, 0: -5: -fa.span[2]),
    yflip=true, color=:isolum, aspect_ratio=2, title="t = $(sol.t[ind])")

snaps = [myplot(i) for i in ind-700: 200: ind+500]

plot(snaps..., layout=(length(snaps), 1), size=(600, 1800))
