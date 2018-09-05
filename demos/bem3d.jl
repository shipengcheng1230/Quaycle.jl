using Distributed
addprocs(4)

@everywhere using SharedArrays

## parameters settings
const ms2mmyr = 365 * 86400 * 1e3
const ρ = 2670.0 # kg/m³
const cs = 3044.0 # m/s
const Vpl = 100.0 # mm/yr
const V0 = 3.2e4 # mm/yr
const f0 = 0.6

# parameters implicit by above
const μ = 0.3 # Bar·km/mm
const λ = μ # poisson material
const α = (λ + μ) / (λ + 2μ)
const η = μ / 2(cs * 1e-3 * 365 * 86400) # Bar·yr/mm

## fault geometry
const dip = 10.0
const depth = 0.0
const lf = 60.0
const nl = 512
const nd = 256

Δl = lf / nl
Δd = Δl
x = collect(range(-lf/2.0+Δl/2, stop=lf/2.0-Δl/2, length=nl))
ξ = collect(range(0., stop=-Δd*(nd-1), step=-Δd)) .- Δd/2.0
z = ξ .* sinpi(dip/180)
y = ξ .* cospi(dip/180)
al = cat(x .- Δl/2, x .+ Δl/2, dims=2)
aw = cat(ξ .- Δd/2, ξ .+ Δd/2, dims=2)
disl = [0., 1., 0.]

## stiffness tensor
@everywhere include(joinpath(dirname(@__DIR__), "src/dc3d.jl"))
const nrept = 2

@everywhere function cal_stiff(x, y, z, α, depth, dip, al, aw, disl, nrept, lf, μ)
    u = zeros(Float64, 12)
    for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i*lf, aw, disl)
    end
    σzz = μ * (3u[12] + u[4] + u[8])
    σyy = μ * (3u[8] + u[4] + u[12])
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180)) / disl[2]
end

function stiff_serial()
    K = zeros(Float64, nl, nd, nd)
    Threads.@threads for l = 1: nd
        for j = 1: nd, i = 1: nl
            K[i,j,l] = cal_stiff(x[i], y[j], z[j], α, depth, dip, al[1,:], aw[l,:], disl, nrept, lf, μ)
        end
    end
    return K
end

function stiff_parfor(nl, nd, x, y, z, α, depth, dip, al, aw, disl, nrept, lf, μ)
    K = SharedArray{Float64}(nl, nd, nd)
    @sync @distributed for l = 1: nd
        for j = 1: nd, i = 1: nl
            K[i,j,l] = cal_stiff(x[i], y[j], z[j], α, depth, dip, al[1,:], aw[l,:], disl, nrept, lf, μ)
        end
    end
end

# @time K = stiff_serial()
@time K = stiff_parfor(nl, nd, x, y, z, α, depth, dip, al, aw, disl, nrept, lf)

## construct a full 2nd order tensor
function full_tensor(K)
    kf = zeros(Float64, nl, nd, nl, nd)
    for t = 1: nd, s = 1: nl, j = 1: nd, i = 1: nl
        kf[i,j,s,t] = K[abs(i-s)+1,j,t]
    end
    kf
end

KF = full_tensor(K)

##
using FFTW
FFTW.set_num_threads(8)
tmp = zeros(Float64, nl)
plan = plan_rfft(tmp, flags=FFTW.MEASURE)

function fft_tensor(K)
    kd = Array{ComplexF64}(undef, Int(nl/2)+1, nd, nd)
    for t = 1: nd, j = 1: nd
        kd[:,j,t] = plan * K[:,j,t]
    end
    return kd
end

KD = fft_tensor(K)

## frictional properties
a = 0.015 .* ones(nl, nd)
b = 0.0115 .* ones(nl, nd)
σ = 150.0 .* ones(nl, nd)
L = 5.0

left_patch = @. -22.5 ≤ x ≤ -2.5
right_patch = @. 2.5 ≤ x ≤ 22.5
vert_patch = @. -20.0 ≤ ξ ≤ -10.0
b[left_patch, vert_patch] .= 0.0185
b[right_patch, vert_patch] .= 0.0185
σ[left_patch, :] .= 15.

const a_ = a
const b_ = b
const σ_ = σ
const L_ = L

## check profile
# using PyPlot
# plt[:imshow](σ)
# plt[:imshow](a-b)

## set ODEs
using DifferentialEquations
using TensorOperations

ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv, relv = [zeros(Float64, nl, nd) for _ in 1: 6]

function f_regularized!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv, relv, ST)
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
    @. relv = Vpl - v
    @tensor begin
        dμ_dt[a,b] = ST[a,b,c,d] * relv[c,d]
    end
    @. dv = (dμ_dt - dμ_dθ * dθ) / (dμ_dv + η)
end

function f_general!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv, relv, ST)
    v = selectdim(u, 3, 1)
    θ = selectdim(u, 3, 2)
    dv = selectdim(du, 3, 1)
    dθ = selectdim(du, 3, 2)

    # make sure θ don't go below 0.0
    clamp!(θ, 0.0, Inf)

    @. dθ = 1 - v * θ / L_
    @. relv = Vpl - v
    # @tensor dμ_dt[a,b] = ST[a,b,c,d] * relv[c,d]
    __dμdt__!(dμ_dt, ST, relv)

    @. dμ_dθ = σ_ * b_ / θ
    @. dμ_dv = σ_ * a_ / v
    @. dv = (dμ_dt - dμ_dθ * dθ) / (dμ_dv + η)
end

function __dμdt__!(out, k, v)
    # @tensor out[a,b] = k[a,b,c,d] * v[c,d]
    fill!(out, 0.)
    @inbounds for l = 1: nd
        for j = 1: nd
            out[:,j] .+= plan \ ((k[:,j,l]) .* (plan * v[:,l]))
        end
    end
end

f1! = (du, u, p, t) -> f_regularized!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv, relv, KF)
f2! = (du, u, p, t) -> f_general!(du, u, p, t, ϕ1, ϕ2, dμ_dt, dμ_dθ, dμ_dv, relv, KD)

v0 = Vpl .* ones(nl, nd)
θ0 = L ./ v0 ./ 1.1
u0 = cat(v0, θ0, dims=3)
tspan = (0.0, 5.0)
prob = ODEProblem(f2!, u0, tspan)
@time sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, progress=true)

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
