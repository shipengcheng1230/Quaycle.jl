##
addprocs(4)

## parameters setting
# pre-defined
const nl = 128 # cell # along strike
const nd = 64 # cell # along dip
const μ = 0.3 # Bar·km/mm
const λ = μ # isotropic medium
const lf = 60.0 # length of frictional law applied [km]
const cs = 3000.0 # m/s
const μ0 = 0.6 # reference frictional coefficient
const nrept = 5 # repetition # of fault length

# implicit
const α = (λ + μ) / (λ + 2μ) # medium constant
const η = μ / 2(cs * 1e-3 * 365 * 86400) # Bar·yr/mm
const Δl = lf / nl # cell size along strike
const Δd = Δl # cell size along dip, here square cells are used
const maxdep = Δd * nd # max depth of where frictional law is applied

## calculate stiffness matrix
dc3d_file = joinpath(dirname(@__DIR__), "src/dc3d.jl")
@everywhere include($dc3d_file)

# fault setting
const dip = 10.0 # dip angle [degree]
const disl = [0., 1., 0.]
const cd = cospi(dip / 180)
const sd = sinpi(dip / 180)
const c2d = cospi(2dip / 180)
const s2d = sinpi(2dip / 180)
const depth = 0.0 # depth of fault origin [km]

# grid coordinates
# centroid of each cell
const x = collect(linspace(0, lf, nl)) - lf / 2
const ξ = collect(linspace(0 - Δd / 2.0, -maxdep, nd))
const y = ξ * cd
const z = ξ * sd
const al = hcat([x - Δl / 2.0 x + Δl / 2.0])
const aw = hcat([ξ + Δd / 2.0 ξ - Δd / 2.0])

KK = SharedArray{Float64}(2nl-1, nd, nd)

# define the kernel
@everywhere function periodic_stiffness(_x, _y, _z, α, depth, dip, _al, _aw, disl, nrept, lenrept, μ)

    @views function net_plane_shear_stress(u)
        σzz = μ * (3u[12] + u[4] + u[8])
        σyy = μ * (3u[8] + u[4] + u[12])
        τyz = μ * (u[9] + u[11])
        tn = (σzz - σyy) * sinpi(2dip / 180) / 2 + τyz * cospi(2dip / 180)
    end

    k = zero(Float64)
    for r = -nrept: nrept
        u = dc3d_okada(_x, _y, _z, α, depth, dip, _al + r * lenrept, _aw, disl)
        tn = net_plane_shear_stress(u)
        k += -tn / disl[2]
    end
    k
end

@everywhere function fill_K_chunk!(K, lrange, x, y, z, α, depth, dip, al, aw, disl, nrept, lenrept, μ)
    @show lrange
    _nl = Int((size(K, 1) + 1) / 2)
    _nd = size(K, 2)

    for l in lrange, j in 1: _nd, i in -_nl + 1: _nl - 1
        if i < 0
            K[i+_nl, j, l] = periodic_stiffness(x[1], y[j], z[j], α, depth, dip, al[1-i, :], aw[l, :], disl, nrept, lenrept, μ)
        else
            K[i+_nl, j, l] = periodic_stiffness(x[end], y[j], z[j], α, depth, dip, al[end-i, :], aw[l, :], disl, nrept, lenrept, μ)
        end
    end
end

@everywhere function local_range(K)
    idx = indexpids(K)
    nchunks = length(procs(K))
    splits = [round(Int, s) for s in linspace(0, size(K, 3), nchunks + 1)]
    splits[idx] + 1: splits[idx+1]
end

@everywhere function fill_K_shared_chunks!(K, x, y, z, α, depth, dip, al, aw, disl, nrept, lenrept, μ)
    fill_K_chunk!(K, local_range(K), x, y, z, α, depth, dip, al, aw, disl, nrept, lenrept, μ)
end

function fill_K_shared!(K, x, y, z, α, depth, dip, al, aw, disl, nrept, lenrept, μ)
    @sync begin
        for p in procs(K)
            @async remotecall_wait(fill_K_shared_chunks!, p, K, x, y, z, α, depth, dip, al, aw, disl, nrept, lenrept, μ)
        end
    end
    K
end

K2 = SharedArray{Float64}(2nl-1, nd, nd)
@time fill_K_shared!(K2, x, y, z, α, depth, dip, al, aw, disl, nrept, lf, μ)



function fill_K!(K)
    @views function net_shear_stress(u)
        σzz = μ * (3u[12] + u[4] + u[8])
        σyy = μ * (3u[8] + u[4] + u[12])
        τyz = μ * (u[9] + u[11])
        tn = (σzz - σyy) * s2d / 2 + τyz * c2d
        kk = -tn / disl[2]
    end

    function periodic_summation(ik, j, l)
        k = zero(Float64)
        for r = -nrept: nrept
            if ik < 0
                u = dc3d_okada(x[1], y[j], z[j], α, depth, dip, al[1-ik, :] + r * lf, aw[l, :], disl)
            else
                u = dc3d_okada(x[end], y[j], z[j], α, depth, dip, al[end-ik, :] + r * lf, aw[l, :], disl)
            end
            k += net_shear_stress(u)
        end
        k
    end

    for l = 1: nd
        for j = 1: nd
            for ik = -nl + 1: nl - 1
                K[ik+nl, j, l] = periodic_summation(ik, j, l)
            end
        end
    end

end

@time fill_K!(KK)
