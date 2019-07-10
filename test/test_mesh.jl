using Test
using JuEQ: triangle_geometric_vector!

@testset "Basic Mesh Generator" begin
    @testset "Line for dc3d" begin
        msh = gen_mesh(Val(:LineOkada), 100.0, 2., 33.0)
        @test msh.nξ == 50
        @test msh.aξ[end][1] == msh.ξ[end] - msh.Δξ/2
        @test msh.y[end] == msh.ξ[end] * cosd(msh.dip)
        @test msh.z[end] == msh.ξ[end] * sind(msh.dip)
    end

    @testset "Rect for dc3d" begin
        msh = gen_mesh(Val(:RectOkada), 100.0, 50.0, 2.0, 2.0, 33.0)
        @test msh.nx == 50
        @test msh.nξ == 25
        @test msh.aξ[end][1] == msh.ξ[end] - msh.Δξ/2
        @test msh.y[end] == msh.ξ[end] * cosd(msh.dip)
        @test msh.z[end] == msh.ξ[end] * sind(msh.dip)
        @test msh.x[end] - msh.x[1] == msh.Δx * (msh.nx - 1)
    end
end

@testset "compute strike and downdip vector" begin
    ss, ds, ts = [zeros(3) for _ in 1: 3]

    @testset "parallel to x" begin
        θ = 43.0
        A = [-10.0, 0.0, 0.0]
        B = [10.0, 0.0, 0.0]
        C = [10.0, -20.0 * cosd(θ), -20.0 * sind(θ)]
        triangle_geometric_vector!(A, B, C, ss, ds, ts)
        @test atand(ds[3] / ds[2]) ≈ θ
        @test atand(ts[2] / ts[3]) ≈ -θ
        @test ss ≈ [1, 0, 0]
    end

    @testset "parallel to y" begin
        θ = 57.0
        A = [0.0, -10.0, 0.0]
        B = [0.0, 10.0, 0.0]
        C = [20.0 * cosd(θ), 10.0, -20.0 * sind(θ)]
        triangle_geometric_vector!(A, B, C, ss, ds, ts)
        @test atand(ds[3] / ds[1]) ≈ -θ
        @test atand(ts[1] / ts[3]) ≈ θ
        @test ss ≈ [0, 1, 0]
    end

    @testset "parallel to z" begin
        A = [10.0, 0.0, 0.0]
        θ = 83.0
        B = [0.0, -A[1] * tand(θ), 0.0]
        C = [B[1], B[2], -20.0]
        triangle_geometric_vector!(A, B, C, ss, ds, ts)
        @test atand(ss[2] / ss[1]) ≈ θ
        @test atand(ts[1] / ts[2]) ≈ -θ
        @test ds == [0, 0, 1]
    end

    @testset "horizontally coplanar" begin
        A = [0.0, -10.0, 0.0]
        B = [10.0, 10.0, 0.0]
        C = [10.0, 0.0, 0.0]
        @test_logs (:warn, "Coplanar at `z = const` where we cannot tell strike and downdip direction.") triangle_geometric_vector!(A, B, C, ss, ds, ts)
        @test ts ≈ [0, 0, 1]
    end

    @testset "swap ordering" begin
        @testset "clockwise" begin
            A = [-10.0, 0.0, 0.0]
            B = [10.0, 0.0, 0.0]
            C = [10.0, -1.0, -5.0]
            @test cross(B-A, C-A)[3] < 0
            triangle_geometric_vector!(A, B, C, ss, ds, ts)
            @test ts[3] > 0
            @test cross(B-A, C-A)[3] > 0
            @test A == [10.0, 0.0, 0.0]
            @test B == [-10.0, 0.0, 0.0]
            @test C == [10.0, -1.0, -5.0]
        end

        @testset "counter clockwise" begin
            A = [-10.0, 0.0, 0.0]
            B = [10.0, 0.0, 0.0]
            C = [10.0, 1.0, -5.0]
            @test cross(B-A, C-A)[3] > 0
            triangle_geometric_vector!(A, B, C, ss, ds, ts)
            @test ts[3] > 0
            @test cross(B-A, C-A)[3] > 0
            @test A == [-10.0, 0.0, 0.0]
            @test B == [10.0, 0.0, 0.0]
            @test C == [10.0, 1.0, -5.0]
        end

        @testset "parallel to z" begin
            A = [-10.0, 0.0, 0.0]
            B = [10.0, 0.0, 0.0]
            C = [10.0, 0.0, -5.0]
            @test cross(B-A, C-A)[2] > 0
            triangle_geometric_vector!(A, B, C, ss, ds, ts)
            @test ts[2] < 0
            @test cross(B-A, C-A)[2] < 0
            @test A == [10.0, 0.0, 0.0]
            @test B == [-10.0, 0.0, 0.0]
            @test C == [10.0, 0.0, -5.0]
        end
    end
end
