using Distributed
addprocs(4)

## parameters settings
const ms2mmyr = 365 * 86400 * 1e3
const ρ = 2670.0 # kg/m³
const cs = 3464.0 # m/s
const Vp = 140 # mm/yr
const V0 = 1e-6 * ms2mmyr # mm/yr
const f0 = 0.6

# parameters implicit by above
const μ = cs^2 * ρ / 1e5 / 1e6 # Bar·km/mm
const λ = μ # poisson material
const α = (λ + μ) / (λ + 2μ)
const η = μ / 2(cs * 1e-3 * 365 * 86400) # Bar·yr/mm

## fault geometry
const dip = 10.0
const depth = 0.0
const lf = 60.0
const nl = 128
const nd = 64

Δl = lf / nl
Δd = Δl
x = collect(range(-lf/2., stop=lf/2., length=nl))
ξ = collect(range(0., stop=-Δd*(nd-1), step=-Δd)) .- Δd/2
z = ξ .* sinpi(dip/180)
y = ξ .* cospi(dip/180)
al = cat(x .- Δl/2, x .+ Δl/2, dims=2)
aw = cat(ξ .- Δd/2, ξ .+ Δd/2, dims=2)
disl = [0., 1., 0.]

## stiffness tensor
include(joinpath(dirname(@__DIR__), "src/dc3d.jl"))
const nrept = 5

function cal_stiff(x, y, z, α, depth, dip, al, aw, disl)
    u = zeros(Float64, 12)
    for i = -nrept: nrept
        u += dc3d_okada(x, y, z, α, depth, dip, al .+ i*lf, aw, disl)
    end
    σzz = μ * (3u[12] + u[4] + u[8])
    σyy = μ * (3u[8] + u[4] + u[12])
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180)) / disl[2]
end

function stiff_serial()
    K = zeros(Float64, nl, nd, nd)
    for l = 1: nd, j = 1: nd, i = 1: nl
        K[i,j,l] = cal_stiff(x[i], y[j], z[j], α, depth, dip, al[1,:], aw[l,:], disl)
    end
    return K
end

@time K = stiff_serial()

## frictional properties
a = 0.015 .* ones(nl, nd)
b = 0.0115 .* ones(nl, nd)
σ = 200.0 .* ones(nl, nd)
L = 5.0

left_patch = @. -22.5 ≤ x < -2.5
right_patch = @. 2.5 < x ≤ 22.5
vert_patch = @. -20.0 ≤ ξ < -10.0
b[left_patch, vert_patch] .= 0.0185
b[right_patch, vert_patch] .= 0.0185
σ[left_patch, :] .= 20.

const a_ = a
const b_ = b
const σ_ = σ
const L_ = L

## check profile
using PyPlot
plt[:imshow](σ)
plt[:imshow](a-b)

## set ODEs
using DifferentialEquations
using ToeplitzMatrices

ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv = [zeros(Float64, nl, nd) for _ in 1: 5]

function f_regularized!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv, stiff)
    v = selectdim(u, 3, 1)
    θ = selectdim(u, 3, 2)
    dv = selectdim(du, 3, 1)
    dθ = selectdim(du, 3, 2)

    # make sure θ don't go below 0.0
    clamp!(θ, 0.0, Inf)

    # using regularized form
    @. ϕ1 = exp((f0 + b_ * log(V0 * θ / L_)) / a_) / 2V0
    @. ϕ2 = σ_ * ϕ1 / sqrt(1 + (v * ϕ1)^2)
    @. dμ_dv = a_ * ϕ2
    @. dμ_dθ = b_ / θ * v * ϕ2
    @. dθ = 1 - v * θ / L_

    fill!(dμ_dt, 0.0)
    for t = 1: nd, s = 1: nl, j = 1: nd
        @simd for i = 1: nl
            @fastmath @inbounds dμ_dt[i,j] += stiff[abs(i-s)+1,j,t] * (Vp - v[s,t])
        end
    end
    @. dv = (dμ_dt - dμ_dθ * dθ) / (dμ_dv + η)
end

f! = (du, u, p, t) -> f_regularized!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv, K)

v0 = Vp .* ones(nl, nd)
θ0 = L ./ v0 ./ 1.1
u0 = cat(v0, θ0, dims=3)
tspan = (0.0, 10.0)
prob = ODEProblem(f!, u0, tspan)
@time sol = solve(prob, OwrenZen5(), reltol=1e-8, abstol=1e-8, progress=true)

## save the results
using FileIO, HDF5
uu = VectorOfArray(sol.u)
u = convert(Array, uu)

v = selectdim(u, 3, 1)
θ = selectdim(u, 3, 2)

_v = v / ms2mmyr
_θ = θ * 365 * 86400

h5open(joinpath(@__DIR__, "bem3d_solution.h5"), "w") do f
    g = g_create(f, "bem3d")
    g["t"] = sol.t # yr
    g["velocity"] = _v # m/s
    g["state"] = _θ # s
end
