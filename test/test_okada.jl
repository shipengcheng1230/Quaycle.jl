using Base.Test
using JuEQ

# checklist from DC3D manual
flag, u, ∇u = dc3d_wrapper([10., 20., -30.], 2./3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
@test flag == 0
@test all(@. isapprox(u, [-37.8981, 63.1789, 14.9607], atol=1e-4))

# positive z given
flag, u, ∇u = dc3d_wrapper([10., 20., 30.], 2./3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
@test flag == 2

# singularity at fault edge
flag, u, ∇u = dc3d_wrapper([-80., 0., -50.], 2./3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
@test flag == 1
