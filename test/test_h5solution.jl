using Test
using HDF5
using Statistics

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
        getu = (u, t, integrator) -> (u.x[1], u.x[2], u.x[3])
        wsolve(prob, Tsit5(), tmp, 50, getu, ["u1", "u2", "u3"], "t")
        @test h5read(tmp, "t") == sol.t
        for m in eachindex(u0.x)
            x = Array(VectorOfArray([sol.u[i].x[m] for i in eachindex(sol.t)]))
            @test h5read(tmp, ustrs[m]) == x
        end
        # @test h5readattr(tmp, "t")["retcode"] == "Success"

        # appended storage
        prob2 = ODEProblem(foo, u0, (6000.0, 7000.0))
        sol2 = solve(prob2, Tsit5())
        wsolve(prob2, Tsit5(), tmp, 50, getu, ["u1", "u2", "u3"], "t"; append=true)
        tcat = cat(sol.t, sol2.t; dims=1)
        ucat = cat(sol.u, sol2.u; dims=1)
        @test tcat == h5read(tmp, "t")
        for m in eachindex(u0.x)
            x = Array(VectorOfArray([ucat[i].x[m] for i in eachindex(tcat)]))
            @test h5read(tmp, ustrs[m]) ≈ x
        end

        # strided storage
        stride = 11
        wsolve(prob, Tsit5(), tmp, 50, getu, ["u1", "u2", "u3"], "t"; stride=stride, force=true)
        @test length(h5read(tmp, "t")) == length(sol.t) ÷ stride + 1
        for m in eachindex(u0.x)
            x = Array(VectorOfArray([sol.u[i].x[m] for i in 1: stride: length(sol.t)]))
            @test h5read(tmp, ustrs[m]) == x
        end

        rm(tmp)
    end
end

@testset "trim solution" begin
    t = collect(1.: 10.)
    v = zeros(5, 10, length(t))
    v[:,:,1:2:end] .= 2
    tmp = tempname()
    tmpout = tempname()
    h5write(tmp, "t", t)
    h5write(tmp, "v", v)
    h5trimsolution(tmp, tmpout, "t", ["v"], u -> mean(u[1]) > 1, t -> t > 2; nstep=10)
    t′ = h5read(tmpout, "t")
    v′ = h5read(tmpout, "v")
    @test length(t′) == 4 # only 3, 5, 7, 9 are qualified
    @test size(v′, 3) == 4
    @test minimum(v′) == 2 # only retain large mean components
    rm(tmp)
    rm(tmpout)
end
