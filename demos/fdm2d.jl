## Erickson et al., 2014

## domain parameters
const Ly = 40e3
const Lz = 40e3
const Δy = 1e3
const Δz = 1e3
const Ny = Int(Ly / Δy)
const Nz = Int(Lz / Δz)

##
using LinearAlgebra

function H_op(dim::Integer)
    m = diagm(0=>ones(Float64, convert(Dims, (dim,))))
    m[1,1] = 0.5
    m[end,end] = 0.5
    return m
end

Hy = Δy .* H_op(Ny+1)
Hz = Δz .* H_op(Nz+1)

function B_op(dim::Integer)
    m = zeros(Float64, dim)
    m[1] = -1.
    m[end] = -1.
    diagm(0=>m)
end

By = B_op(Ny+1)
Bz = B_op(Nz+1)

function S_op(dim::Integer)
    m = Matrix{Float64}(I, dim, dim)
    m[1,1:3] = [-1.5, 2, -0.5]
    m[end,end-2:end] = [0.5, -2, 1.5]
    return m
end

Sy = (1/Δy) .* S_op(Ny+1)
Sz = (1/Δz) .* S_op(Nz+1)

Iy = Matrix{Float64}(I, Ny+1, Ny+1)
Iz = Matrix{Float64}(I, Nz+1, Nz+1)

function D_op(dim::Integer)
    diag = 0.5 .* ones(Float64, dim-1)
    m = diagm(1=>diag, -1=>-diag)
    m[1,1] = -1.
    m[end,end] = 1.
    m[end,end-1] = -1.
    m[1,2] = 1.
    return m
end

Dy = (1/Δy) .* D_op(Ny+1)
Dz = (1/Δz) .* D_op(Nz+1)

function D2_op(dim::Integer)
    m = diagm(0=>-2 .*ones(Float64, dim), 1=>ones(Float64, dim-1), -1=>ones(Float64, dim-1))
    m[1,1:3] = [1., -2., 1.]
    m[end,end-2:end] = [1., -2., 1.]
    return m
end

D2y = (1/Δy^2) .* D2_op(Ny+1)
D2z = (1/Δz^2) .* D2_op(Nz+1)

function R_op()

end

function C2_op(dim::Integer)
    m = Matrix{Float64}(I, dim, dim)
    m[1,1] = 0.
    m[end,end] = 0.
    return m
end

C2y = C2_op(Ny+1)
C2z = C2_op(Nz+1)

Rμy = (0.25 * Δy^7) .* (kron(D2y', Iz)) * (C)
