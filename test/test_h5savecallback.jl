using Test
using HDF5

@testset "h5savecallback" begin
    @testset "vector output" begin
        function lorenz(du, u, p, t)
            du[1] = 10.0 * (u[2] - u[1])
            du[2] = u[1] * (28.0 - u[3]) - u[2]
            du[3] = u[1] * u[2] - (8/3) * u[3]
        end

        tmp = tempname()
        cbfun = @h5savecallback(tmp, 100.0, 10, (3,), Float64)
        u0 = [1.0; 0.0; 0.0]
        tspan = (0.0, 100.0)
        prob = ODEProblem(lorenz, u0, tspan)
        sol = solve(prob, Tsit5(); save_everystep=false, callback=cbfun)
        sol2 = solve(prob, Tsit5())
        res_u = h5read(tmp, "u")
        u1s = res_u[2, :]
        u2s = [x[2] for x in sol2.u]
        all(u1s .== u2s)
        t1 = h5read(tmp, "t")
        t2 = sol2.t
        @test sum(t1 .- t2) < 1e-8
        rm(tmp)
    end

    @testset "matrix output" begin
        A  = [1. 0  0 -5
              4 -2  4 -3
             -4  0  0  1
              5 -2  2  3]
        tmp = tempname()
        cbfun = @h5savecallback(tmp, 100.0, 25, (4, 2), Float64)
        u0 = rand(4, 2)
        tspan = (0.0,100.0)
        f(u, p, t) = A * u
        prob = ODEProblem(f, u0, tspan)
        sol = solve(prob, Tsit5(); save_everystep=false, callback=cbfun)
        sol2 = solve(prob, Tsit5())
        res_u = h5read(tmp, "u")
        u1s = res_u[2,1,:]
        u2s = [x[2,1] for x in sol2.u]
        @test u1s == u2s
        t1 = h5read(tmp, "t")
        t2 = sol2.t
        @test sum(t1 .- t2) < 1e-8
        rm(tmp)
    end

    @testset "ArrayPartition" begin
        function foo(du, u, p, t)
            du.x[1] .= u.x[1] / 10
            du.x[2] .= u.x[2] / 100
            du.x[3] .= u.x[3] / 1000
        end
        u0 = ArrayPartition(rand(2, 3), rand(5), rand(3, 2))
        tspan = (0.0, 5000.0)
        prob = ODEProblem(foo, u0, tspan)
        sol = solve(prob, Tsit5())
        ustrs = ["u1", "u2", "u3"]
        tmp = tempname()
        B = create_h5buffer(tmp, u0, 100, tspan[2], ustrs, "t")
        getu = (u, t, integrator, b) -> (u.x[1], u.x[2], u.x[3])
        callback = (u, t, integrator) -> h5savebuffercallback_kernel(u, t, integrator, B, getu)
        fcb = FunctionCallingCallback(callback)
        solve(prob, Tsit5(); callback=fcb)
        @test h5read(tmp, "t") == sol.t
        for m in 1: length(u0.x)
            x = Array(VectorOfArray([sol.u[i].x[m] for i in 1: length(sol.t)]))
            @test h5read(tmp, ustrs[m]) == x
        end
        rm(tmp)
    end
end
