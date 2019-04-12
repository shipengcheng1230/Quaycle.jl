# !!! note
#     This example is adapted from [Wei, 2016 AGU](https://agu.confex.com/agu/fm16/meetingapp.cgi/Paper/132124)

# !!! tip
#     It will automatically use parallel scheme if `nprocs() ≂̸ 1` when building stiffness tensor. To do so:
#     ```julia linenums="1"
#     using Distributed
#     addprocs(4); # add # of cores you desire
#     using JuEQ
#     ```

# ### Define Parameters
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

fa = fault(Val(:CSFS), STRIKING(), (80., 10.), (0.5, 0.5))

# Generate grids:


# !!! tip
#     It is recommended (from Yajing Liu's personal communication) to add buffer zones adjacent the horizontal edges
#     to immitate *zero* dislocation at the ridge region.
#     Basically, it affects how the stiffness tensor are periodically summed. To what extent it alters the results remains further testing.

#     Under the hood, it shall impose buffer areas on both sides of along-strike, each of which has a length of `bufferratio/2*fa[:x]`.
#     Thus, the stiffness contributions falling into those buffer zone shall be neglected, which is equivalent to impose zero-slip correspondingly.

# Time for us to establish frictional parameters profile:

frprop = init_friction_prop(fa)
faprop = init_fault_prop(λ, μ, η, vpl, f0, v0)

fill!(frprop.a, 0.015)
fill!(frprop.b, 0.0115)
fill!(frprop.L, 12.0)

left_patch = @. -25. ≤ fa.mesh.x ≤ -5.
right_patch = @. 5. ≤ fa.mesh.x ≤ 25.
vert_patch = @. -6. ≤ fa.mesh.ξ ≤ -1

frprop.b[xor.(left_patch, right_patch), vert_patch] .= 0.0185

σmax = 500.
σ = [min(σmax, 15. + 180. * z) for z in -fa.mesh.ξ]
σ = Matrix(repeat(σ, 1, fa.mesh.nx)')
L = 12.;

frprop.σ .= σ

frprop
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

# ### Construct Model
# Construct our material property profile:

mp = properties(fa, gd, [:a=>a, :b=>b, :L=>L, :σ=>σ, :η=>η, :k=>[:λ=>λ, :μ=>μ], :vpl=>vpl, :f0=>f0, :v0=>v0]);

# Provide the initial condition:

vinit = vpl .* ones(fa.mesh.nx, fa.mesh.nξ)
θ0 = L ./ vinit ./ 1.1
u0 = cat(vinit, θ0, zeros(Float64, fa.mesh.nx, fa.mesh.nξ), dims=3);

prob = assemble(Val(:okada), fa, faprop, frprop, u0, (0., 18.))

# Get our ODEs problem:

# ### Solve Model
# Solve the model:

@time sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6);
sol.t[end]
f = prob.f
prob.p
du = similar(u0)
@time f(du, u0, prob.p, 1.0)

# ### Sanity Check of Results
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
