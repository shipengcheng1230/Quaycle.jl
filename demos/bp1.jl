# This is a demo for SEAS BP1, inspired by Yajing Liu's code
# 2D vertical anti-shear-plain

## let's do parallel
addprocs(4)
# use dc3d
dc3d_file = joinpath(dirname(@__DIR__), "src/dc3d.jl")
# need string insertion here to work around
@everywhere include($dc3d_file)

## parameters are set using BP1 requirements
const ρ = 2670.0
const g = 9.8
const cs = 3464.0
const σ = 50.0e6
const a0 = 0.010
const amax = 0.025
const b0 = 0.015
const Dc = 0.008
const Vp = 1e-9
const Vinit = 1e-9
const V0 = 1e-6
const f0 = 0.6
const H = 15.0e3
const h = 3.0e3
const Wf = 40.0e3
const Δz = 25.0
const tf = 3000.0

# parameters implicit by above
const μ = cs^2 * ρ
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

fill_stiff_matrix!(K)
## modeling
