using Test
using DifferentialEquations

@testset "h5save" begin
    @testset "vector output" begin
        function lorenz(du, u, p, t)
            du[1] = 10.0 * (u[2] - u[1])
            du[2] = u[1] * (28.0 - u[3]) - u[2]
            du[3] = u[1] * u[2] - (8/3) * u[3]
        end

        tmp = tempname()
        cbfun = @h5save(tmp, 100.0, 10, (3,), 1, Float64)

        u0 = [1.0; 0.0; 0.0]
        tspan = (0.0, 100.0)
        prob = ODEProblem(lorenz, u0, tspan)
        sol = solve(prob; save_everystep=false, callback=FunctionCallingCallback(cbfun; func_everystep=true), funcat=[100.0])
        sol2 = solve(prob)
        res_u = h5read(tmp, "u")
        u1s = res_u[2, :]
        u2s = [x[2] for x in sol2.u]
        all(u1s .== u2s)
        t1 = h5read(tmp, "t")
        t2 = sol2.t
        @test sum(t1 .- t2) < 1e-8
    end

    @testset "matrix output" begin
        A  = [1. 0  0 -5
              4 -2  4 -3
             -4  0  0  1
              5 -2  2  3]
        tmp = tempname()
        cbfun = @h5save(tmp, 100.0, 25, (4, 2), 2, Float64)
        u0 = rand(4, 2)
        tspan = (0.0,100.0)
        f(u, p, t) = A * u
        prob = ODEProblem(f, u0, tspan)
        sol = solve(prob; save_everystep=false, callback=FunctionCallingCallback(cbfun; func_everystep=true), funcat=[100.0])
        sol2 = solve(prob)
        res_u = h5read(tmp, "u")
        u1s = res_u[2,1,:]
        u2s = [x[2,1] for x in sol2.u]
        @test u1s == u2s
        t1 = h5read(tmp, "t")
        t2 = sol2.t
        @test sum(t1 .- t2) < 1e-8
    end
end
