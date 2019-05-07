# greens function ϵ vs u for tet4 element
export sbarbot_disp_tet4, sbarbot_disp_tet4!

function sbarbot_disp_tet4(quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, nu::R
    ) where {R, U, Q}

    u = Vector{R}(undef, 3)
    sbarbot_disp_tet4!(u, quadrature, x1, x2, x3, A, B, C, D, e11, e12, e13, e22, e23, e33, nu)
    return u
end

"""
                      / North (x1)
                     /
        surface     /
      -------------+-------------- East (x2)
                  /|
                 / |     + A
                /  |    /  .
                   |   /     .
                   |  /        .
                   | /           .
                   |/              + B
                   /            .  |
                  /|          /    |
                 / :       .       |
                /  |    /          |
               /   : .             |
              /   /|               |
             / .   :               |
            +------|---------------+
          C        :                 D
"""
function sbarbot_disp_tet4!(u::W, quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, nu::R
    ) where {R, U, Q, W}

    lambda = 2 * nu / (1 - 2 * nu)
    ekk = e11 + e22 + e33

    nA = cross(C - B, D - B)
    nB = cross(D - C, A - C)
    nC = cross(A - D, B - D)
    nD = cross(B - A, C - A)

    nA /= norm(nA)
    nB /= norm(nB)
    nC /= norm(nC)
    nD /= norm(nD)

    if nA' * (A .- (B .+ C .+ D) ./ 3) > zero(R); nA *= -one(R) end
    if nB' * (B .- (A .+ C .+ D) ./ 3) > zero(R); nB *= -one(R) end
    if nC' * (C .- (B .+ A .+ D) ./ 3) > zero(R); nC *= -one(R) end
    if nD' * (D .- (B .+ C .+ A) ./ 3) > zero(R); nD *= -one(R) end

    ABC = norm(cross(C .- A, B .- A)) / 2
    BCD = norm(cross(D .- B, C .- B)) / 2
    CDA = norm(cross(A .- C, D .- C)) / 2
    DAB = norm(cross(B .- D, A .- D)) / 2

    m11 = lambda * ekk + 2 * e11
    m12 = 2 * e12
    m13 = 2 * e13
    m22 = lambda * ekk + 2 * e22
    m23 = 2 * e23
    m33 = lambda * ekk + 2 * e33

    let lambda=lambda, x1=x1, x2=x2, x3=x3, A=A, B=B, C=C, D=D, e11=e11, e12=e12, e13=e13, e22=e22, e23=e23, e33=e33, nu=nu, nA=nA, nB=nB, nC=nC, ABC=ABC, BCD=BCD, CDA=CDA, DAB=DAB, m11=m11, m12=m12, m13=m13, m22=m22, m23=m23, m33=m33

        function r1(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 - y3)
        end

        function r2(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 + y3)
        end

        function G11(y1::R, y2::R, y3::R) where R
            (1 / (16 * pi * (1 - nu)) * (
            (3 - 4 * nu) / r1(y1, y2, y3) + 1 / r2(y1, y2, y3) + (x1 - y1) ^ 2 / r1(y1, y2, y3) ^ 3
            + (3 - 4 * nu) * (x1 - y1) ^ 2 / r2(y1, y2, y3) ^ 3 + 2 * x3 * y3 * (r2(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 5
            + 4 * (1 - 2 * nu) * (1 - nu) * (r2(y1, y2, y3) ^ 2 - (x1 - y1) ^ 2 + r2(y1, y2, y3) * (x3 + y3)) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G12(y1::R, y2::R, y3::R) where R
            ((x1 - y1) * (x2 - y2) / (16 * pi * (1 - nu)) * (
            1  / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) / r2(y1, y2, y3) ^ 3 - 6 * x3 * y3 / r2(y1, y2, y3) ^ 5
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G13(y1::R, y2::R, y3::R) where R
            ((x1 - y1) /(16*pi*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            -6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 + 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G21(y1::R, y2::R, y3::R) where R
            ((x1 - y1) * (x2 - y2) / (16 * pi * (1 - nu)) * (
            1  / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) / r2(y1, y2, y3) ^ 3 - 6 * x3 * y3 / r2(y1, y2, y3) ^ 5
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G22(y1::R, y2::R, y3::R) where R
            (1 / (16 * pi * (1 - nu)) * (
            (3 - 4 * nu) / r1(y1, y2, y3) + 1  / r2(y1, y2, y3) + (x2 - y2) ^ 2  / r1(y1, y2, y3) ^ 3
            +(3 - 4 * nu) * (x2 - y2) ^ 2  / r2(y1, y2, y3) ^ 3 + 2 * x3 * y3 * (r2(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 5
            +4 * (1 - 2 * nu) * (1 - nu) * (r2(y1, y2, y3) ^ 2 - (x2 - y2) ^ 2 + r2(y1, y2, y3) * (x3 + y3)) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G23(y1::R, y2::R, y3::R) where R
            ((x2 - y2) /(16*pi*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            -6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 + 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G31(y1::R, y2::R, y3::R) where R
            ((x1 - y1) /(16*pi*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            +6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 - 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G32(y1::R, y2::R, y3::R) where R
            ((x2 - y2) /(16*pi*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            +6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 - 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G33(y1::R, y2::R, y3::R) where R
            (1 / (16 * pi * (1 - nu)) * (
            (3 - 4 * nu) / r1(y1, y2, y3) + (5 - 12 * nu + 8 * nu ^ 2) / r2(y1, y2, y3) + (x3 - y3) ^ 2  / r1(y1, y2, y3) ^ 3
            +6 * x3 * y3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 5 + ((3 - 4 * nu) * (x3 + y3) ^ 2 - 2 * x3 * y3) / r2(y1, y2, y3) ^ 3))
        end

        function y(u::R, v::R, A::R, B::R, C::R) where R
            A * (1 - u) * (1 - v) / 4 + B * (1 + u) * (1 - v) / 4 + C * (1 + v) / 2
        end

        function IU1(u::R, v::R) where R
            (ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G11(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G21(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G31(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G11(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G21(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G31(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G11(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G21(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G31(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G11(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G21(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G31(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU2(u::R, v::R) where R
            (ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G12(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G22(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G32(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G12(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G22(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G32(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G12(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G22(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G32(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G12(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G22(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G32(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU3(u::R, v::R) where R
            (ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G13(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G23(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G33(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G13(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G23(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G33(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G13(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G23(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G33(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G13(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G23(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G33(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        fill!(u, zero(R))
        x, w = quadrature
        @inbounds @fastmath for k = 1: length(x)
            @simd for j = 1: length(x)
                u[1] += w[j] * w[k] * (1 - x[k]) * IU1(x[j], x[k])
                u[2] += w[j] * w[k] * (1 - x[k]) * IU2(x[j], x[k])
                u[3] += w[j] * w[k] * (1 - x[k]) * IU3(x[j], x[k])
            end
        end

        return nothing
    end

end
