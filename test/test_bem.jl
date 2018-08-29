using Test

@testset "Fault uniform FD grids" begin
    fa = fault(ThrustFault, 10., (60., 30.))

    @testset "Interval corretion" begin
        gd = discretize(fa, fa[:x]/128, fa[:ξ]/64)
        @test abs(gd.x[2] - gd.x[1]) ≈ fa[:x] / 128
        @test abs(gd.ξ[2] - gd.ξ[1]) ≈ fa[:ξ] / 64
    end
end

using Distributed
using .BEM
addprocs(4)
nprocs()
using Profile
using InteractiveUtils

fa = fault(StrikeSlipFault, 40.)
gd = discretize(fa, fa[:ξ]/1600)
ep = HomogeneousElasticProperties(λ=0.32, μ=0.32, cs=3.464*365*86400)
@time stiffness_tensor(fa, gd, ep, ncpus=1)

shear_traction(StrikeSlipFault, ones(12), 1., 1., 1.)

##
x, y, z, dep = 0., gd.ξ .* cospi(fa.dip / 180), gd.ξ .* sinpi(fa.dip / 180), 0.
disl = [1., 0., 0.]
ax = [-2000., 2000.]
