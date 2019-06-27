# !!! note
#     This example is adapted from [Wei, 2016 AGU](https://agu.confex.com/agu/fm16/meetingapp.cgi/Paper/132124)

# !!! tip
#     It will automatically use parallel scheme if `nprocs() ≂̸ 1` when building stiffness tensor. To do so:
#     ```julia linenums="1"
#     using Distributed
#     addprocs(4); # add # of cores you desire
#     @everywhere using JuEQ
#     ```

# First, list all the essential parameters:

using JuEQ
using Plots

ms2mmyr = 365 * 86400 * 1e3
ρ = 2670.0 # kg/m³
cs = 3044.0 # m/s
vpl = 100.0 # mm/yr
v0 = 3.2e4 # mm/yr
f0 = 0.6;
μ = 0.3 # Bar·km/mm
λ = μ # poisson material
α = (λ + μ) / (λ + 2μ)
η = μ / 2(cs * 1e-3 * 365 * 86400); nothing # Bar·yr/mm

# First, create a fault mesh, specify fault type and compute the Green's function.

ft = STRIKING()
mesh = gen_mesh(Val(:RectOkada), 80., 10., 0.5, 0.5, 90.0)
gf = okada_stress_gf_tensor(mesh, λ, μ, ft; buffer_ratio=1); nothing

# Next, establish frictional and fault space parameters:
a = ones(mesh.nx, mesh.nξ) .* 0.015
b = ones(mesh.nx, mesh.nξ) .* 0.0115
L = ones(mesh.nx, mesh.nξ) .* 12.0

left_patch = @. -25. ≤ mesh.x ≤ -5.
right_patch = @. 5. ≤ mesh.x ≤ 25.
vert_patch = @. -6. ≤ mesh.z ≤ -1

b[xor.(left_patch, right_patch), vert_patch] .= 0.0185

σmax = 500.
σ = [min(σmax, 15. + 180. * z) for z in -mesh.z]
σ = Matrix(repeat(σ, 1, mesh.nx)')
prop = RateStateQuasiDynamicProperty(a=a, b=b, L=L, σ=σ, vpl=vpl, f0=f0, v0=v0, η=η); nothing

# Make sure our profile match our expectation:

p1 = plot((a .- b)', seriestype=:heatmap,
    xticks=(collect(1: 40: mesh.nx+1), [-40, -20, 0, 20, 40]),
    yticks=(collect(1: 5: mesh.nξ+1), [0, 5, 10, 15, 20]),
    yflip=true, color=:isolum, aspect_ratio=2, title="a-b",
    );

p2 = heatmap(σ',
    xticks=(collect(1: 40: mesh.nx+1), [-40, -20, 0, 20, 40]),
    yticks=(collect(1: 5: mesh.nξ+1), [0, 5, 10, 15, 20]),
    yflip=true, color=:isolum, aspect_ratio=2, title="\\sigma"
    );

plot(p1, p2, layout=(2, 1))

# Then, provide the initial condition and assemble the ODEs:

vinit = vpl .* ones(mesh.nx, mesh.nξ)
θ0 = L ./ vinit ./ 1.1
u0 = ArrayPartition(vinit, θ0)
prob = assemble(mesh, gf, prop, u0, (0., 18.)); nothing

# !!! tip
#     It is recommended (from Yajing Liu's personal communication) to add buffer zones adjacent the horizontal edges
#     to immitate *zero* dislocation at the ridge region.
#     Basically, it affects how the stiffness tensor are periodically summed. To what extent it alters the results remains further testing.

#     Under the hood, it shall impose buffer areas on both sides of along-strike, each of which has a length of `bufferratio/2*fa[:x]`.
#     Thus, the stiffness contributions falling into those buffer zone shall be neglected, which is equivalent to impose zero-slip correspondingly.

# Afterwards, solve ODEs problem:

sol = solve(prob, VCABM5(), reltol=1e-5, abstol=1e-3); nothing

# Last, take a look at the max velocity time series:

maxv = JuEQ.max_velocity(sol)
plot(sol.t, log10.(maxv / ms2mmyr), xlabel="Time (year)", ylabel="Max Velocity (log10 (m/s))", label="")

# And view some snapshots of ruptures (quasi-dynamic) patterns:

ind = argmax(maxv)
myplot = (ind) -> heatmap(log10.(sol.u[ind].x[1]./ms2mmyr)',
    xticks=(collect(1: 40: mesh.nx+1), [-40, -20, 0, 20, 40]),
    yticks=(collect(1: 5: mesh.nξ+1), [0, 5, 10, 15, 20]),
    yflip=true, color=:isolum, aspect_ratio=2, title="t = $(sol.t[ind])")
snaps = myplot(ind+300)
plot(snaps)
