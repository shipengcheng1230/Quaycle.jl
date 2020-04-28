begin
    using Revise
    # using Quaycle
    using GmshTools
    using Gmsh_SDK_jll
    using BenchmarkTools
    using LinearAlgebra

    Threads.nthreads()
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 40000
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10.0
    ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
end

# begin
#     mesh = gen_mesh(Val(:RectOkada), 10., 10., 10.0/512, 10.0/64, 90.)
#     # gf = stress_greens_func(mesh, 1.0, 1.0, STRIKING(); buffer_ratio=1.0)
#     gf = rand(ComplexF64, mesh.nx, mesh.nξ, mesh.nξ)
#     p = RateStateQuasiDynamicProperty([rand(mesh.nx, mesh.nξ) for _ in 1: 4]..., rand(4)...)
#     u0 = ArrayPartition([rand(mesh.nx, mesh.nξ) for _ in 1: 3]...)
#     prob = assemble(gf, p, u0, (0., 1.0))
#     du = similar(u0)
# end
#
# @benchmark prob.f(du, u0, prob.p, 1.0) samples=40000

begin
    mf = gen_mesh(Val(:RectOkada), 10., 10., 10.0/512, 10.0/64, 90.)
    nx = ny = nz = 10
    nume = prod([nx, ny, nz])
    ϵcomp = (:xx, :xy, :xz, :yy, :yz, :zz)
    σcomp = (:xx, :xy, :xz, :yy, :yz, :zz)
    ee = rand(ComplexF64, mf.nx, mf.nξ, mf.nξ)
    ev = rand(nx * ny * nz * 6, mf.nx * mf.nξ)
    ve = Array(ev') |> copy
    vv = rand(nx * ny * nz * 6, nx * ny * nz * 6)
    gg = Quaycle.ViscoelasticCompositeGreensFunction(ee, ev, ve, vv, nx * ny * nz, 6, 6, ϵcomp, σcomp)

    pe = RateStateQuasiDynamicProperty([rand(mf.nx, mf.nξ) for _ in 1: 4]..., rand(4)...)
    pdisl = DislocationCreepProperty([rand(nume) for _ in 1: 10]...)
    pdiff = DiffusionCreepProperty([rand(nume) for _ in 1: 11]...)
    pc = compose(pe, rand(6), pdisl, pdiff)
    v0 = rand(mf.nx, mf.nξ)
    θ0 = rand(mf.nx, mf.nξ)
    ϵ0 = rand(nume, 6)
    σ0 = rand(nume, 6)
    δ0 = rand(mf.nx, mf.nξ)
    u0 = ArrayPartition(v0, θ0, ϵ0, σ0, δ0)
    prob, alloc = assemble(gg, pc, u0, (0., 1.0); return_alloc=true)
    du = similar(u0)
end

@benchmark prob.f(du, u0, prob.p, 1.0) samples=40000
