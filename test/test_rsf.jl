using Test
using DifferentialEquations

mp = ZeroDimMaterialProfile(ev=DieterichStateLaw, fform=CForm, a=0.001, b=0.0015, L=1e-4, k=50.0, vpl=1e-5, f0=0.6, v0=1e-6, η=0.0, σ=1.0)
prob = EarthquakeCycleProblem(mp, [1e-6, mp.L/1e-6], (0., 100.,))
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
μ = friction.(mp, sol.u)
