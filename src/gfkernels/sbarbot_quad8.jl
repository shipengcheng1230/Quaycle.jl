# greens function ϵ vs u for quad8 element
export sbarbot_disp_quad8, sbarbot_disp_quad8!
export sbarbot_strain_quad8, sbarbot_strain_quad8!
export sbarbot_stress_quad8, sbarbot_stress_quad8!

function sbarbot_disp_quad8(
    x1::R, x2::R, x3::R, q1::R, q2::R, q3::R,
    L::R, T::R, W::R, theta::R,
    epsv11p::R, epsv12p::R, epsv13p::R, epsv22p::R, epsv23p::R, epsv33p::R,
    G::R, nu::R,
    ) where R

    u = Vector{R}(undef, 3)
    sbarbot_disp_quad8!(u, x1, x2, x3, q1, q2, q3, L, T, W, theta, epsv11p, epsv12p, epsv13p, epsv22p, epsv23p, epsv33p, G, nu)
    return u
end

function sbarbot_strain_quad8(
    x1::R, x2::R, x3::R, q1::R, q2::R, q3::R,
    L::R, T::R, W::R, theta::R,
    epsv11p::R, epsv12p::R, epsv13p::R, epsv22p::R, epsv23p::R, epsv33p::R,
    G::R, nu::R,
    ) where R

    ϵ = Vector{R}(undef, 6)
    sbarbot_strain_quad8!(ϵ, x1, x2, x3, q1, q2, q3, L, T, W, theta, epsv11p, epsv12p, epsv13p, epsv22p, epsv23p, epsv33p, G, nu)
    return ϵ
end

function sbarbot_stress_quad8!(σ::AbstractVector,
    x1::R, x2::R, x3::R, q1::R, q2::R, q3::R,
    L::R, T::R, W::R, theta::R,
    epsv11p::R, epsv12p::R, epsv13p::R, epsv22p::R, epsv23p::R, epsv33p::R,
    G::R, nu::R,
    ) where R

    lambda = G * 2 * nu / (1 - 2 * nu)
    sbarbot_strain_quad8!(σ, x1, x2, x3, q1, q2, q3, L, T, W, theta, epsv11p, epsv12p, epsv13p, epsv22p, epsv23p, epsv33p, G, nu)
    ekk = σ[1] + σ[4] + σ[6]
    σ[1] = lambda * ekk + 2G * σ[1]
    σ[2] *= 2G
    σ[3] *= 2G
    σ[4] = lambda * ekk + 2G * σ[4]
    σ[5] *= 2G
    σ[6] = lambda * ekk + 2G * σ[6]
    return nothing
end

function sbarbot_stress_quad8(
    x1::R, x2::R, x3::R, q1::R, q2::R, q3::R,
    L::R, T::R, W::R, theta::R,
    epsv11p::R, epsv12p::R, epsv13p::R, epsv22p::R, epsv23p::R, epsv33p::R,
    G::R, nu::R,
    ) where R

    σ = Vector{R}(undef, 6)
    sbarbot_stress_quad8!(σ, x1, x2, x3, q1, q2, q3, L, T, W, theta, epsv11p, epsv12p, epsv13p, epsv22p, epsv23p, epsv33p, G, nu)
    return σ
end

@doc raw"""
## copyright:
    https://bitbucket.org/sbarbot

                      N (x1)
                     /
                    /| strike (theta)          E (x2)
        q1,q2,q3 ->@--------------------------+
                   |                        w |     +
                   |                        i |    /
                   |                        d |   / s
                   |                        t |  / s
                   |                        h | / e
                   |                          |/ n
                   +--------------------------+  k
                   :       l e n g t h       /  c
                   |                        /  i
                   :                       /  h
                   |                      /  t
                   :                     /
                   |                    +
                   Z (x3)

"""
function sbarbot_disp_quad8!(u::AbstractVector{<:R},
    x1::R, x2::R, x3::R, q1::R, q2::R, q3::R,
    L::R, T::R, W::R, theta::R,
    epsv11p::R, epsv12p::R, epsv13p::R, epsv22p::R, epsv23p::R, epsv33p::R,
    G::R, nu::R,
    ) where R

    lambda = G * 2 * nu / (1 - 2 * nu)
    epsvkk = epsv11p + epsv22p + epsv33p

    t1 = (x1 - q1) * cosd(theta) + (x2 - q2) * sind(theta)
    x2 = -(x1 - q1) * sind(theta) + (x2 - q2) * cosd(theta)
    x1 = t1

    # https://github.com/JuliaLang/julia/issues/15276
    let lambda=lambda, epsvkk=epsvkk, x1=x1, x2=x2, x3=x3, q1=q1, q2=q2, q3=q3, L=L, T=T, W=W, theta=theta, epsv11p=epsv11p, epsv12p=epsv12p, epsv13p=epsv13p, epsv22p=epsv22p, epsv23p=epsv23p, epsv33p=epsv33p, G=G, nu=nu

        function r1(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 - y3)
        end

        function r2(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 + y3)
        end

        function J1112(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
            -1) -4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)^ (
              -1)* (x2 - y2)) -x3* atan(x3, x1 - y1) - 3 * x3*
            atan(3 * x3, x1 - y1) + 4 * nu* x3* atan(-nu* x3, x1 -
              y1) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan(r2(y1, y2, y3)* (-x1 + y1), (
                x2-y2)* (x3 + y3)) -4 * ((-1) + nu)* (x3 - y3)* atan(r1(y1, y2, y3)* (
                  x3 - y3), (x1 - y1)* (x2 - y2)) +3 * y3* atan((-3)* y3,
                    x1 - y1) - y3* atan(y3, x1 - y1) - 4 * nu* y3* atan(
                      nu* y3, x1 - y1) - 4 * ((-1) + nu)* (x3 + y3)* atan(r2(y1, y2, y3)* (x3 + y3), (
                        x1-y1)* (x2 - y2)) +xlogy(-((-3) + 4 * nu)* (x1 -
                          y1), r1(y1, y2, y3) + x2 - y2) +xlogy((5 + 4 * nu* ((-3) + 2 * nu))* (x1 - y1),
                            r2(y1, y2, y3) + x2 - y2) + xlogy((-4)* ((-1) + nu)* (x2 - y2), r1(y1, y2, y3) + x1 -
                              y1) + xlogy((-4)* ((-1) + nu)* (x2 - y2), r2(y1, y2, y3) + x1 - y1)))
        end

        function J1113(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* (x1 + (
            -1)* y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (-((-1) +
            nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
            y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
            ) +x2* atan(-x2, x1 - y1) - 3 * x2* atan(3 * x2, x1 -
              y1) + 4 * nu* x2* atan(-nu* x2, x1 - y1) - 4 * ((-1) + nu)* (
            x2 - y2)* atan(r1(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 - y3)
            ) + 4 * ((-1) + nu)* (x2 - y2)* atan(r2(y1, y2, y3)* (x2 - y2), (x1 -
              y1)* (x3 + y3)) +3 * y2* atan((-3)* y2, x1 - y1) - y2* atan(
                y2, x1 - y1) - 4 * nu* y2* atan(nu* y2, x1 - y1) + xlogy((-1)
            * ((-3) + 4 * nu)* (x1 - y1), r1(y1, y2, y3) + x3 - y3) + xlogy(-(3
                  -6 * nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +xlogy((-4)* ((-1) + nu)* (
                    x3 - y3), r1(y1, y2, y3) + x1 - y1) +xlogy(4 * ((-1) + nu)* (x3 + y3), r2(y1, y2, y3) + x1 + (
                      -1)* y1)))
        end

        function J1123(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)* ((
              x1- y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* ((x1 + (-1)
            * y1)^ 2 + (x3 + y3)^ 2)^ (-1)* (x3* ((x3^ 2 + (x1 - y1)^ 2)* (
            x3^ 2 + (x1 - y1)^ 2 + (x2 - y2)^ 2) +x3* (3 * x3^ 2 + 2 * (x1 + (-1)
            * y1)^ 2 + (x2 - y2)^ 2)* y3 + 3 * x3^ 2 * y3^ 2 + x3* y3^ 3) -((
            - 1) +nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)
            ^ 2) +((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)* y3* (2 * x3 + y3)* ((x1 - y1)
            ^ 2 + (x3 + y3)^ 2)) +2 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)* atan((
              x1-y1)* (x2 - y2)^ (-1)) +x1* atan(-x1, x2 - y2)
            -3 * x1* atan(3 * x1, x2 - y2) + 4 * nu* x1* atan(-nu* x1,
              x2 - y2) + 3 * y1* atan((-3)* y1, x2 - y2) - y1* atan(
                y1, x2 - y2) - 4 * nu* y1* atan(nu* y1, x2 - y2) + 2 * ((-1) +
            2 * nu)* (x1 - y1)* atan(r1(y1, y2, y3)* (-x1 + y1), (x2 - y2)* (x3 +
              (-1)* y3)) +2 * (1 - 2 * nu)^ 2 * (x1 - y1)* atan(r2(y1, y2, y3)* (-x1 +
                y1), (x2 - y2)* (x3 + y3)) +xlogy((-2)* x3, r2(y1, y2, y3) - x2 + y2) + xlogy((
            -1)* ((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x3 - y3) +xlogy(-(3 + (
                  -6)* nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +xlogy(-((-3) + 4 *
                    nu)* (x3 - y3), r1(y1, y2, y3) + x2 - y2) +xlogy(-(5 + 4 * nu* ((-3) + 2 *
                      nu))* (x3 + y3), r2(y1, y2, y3) + x2 - y2)))
        end

        function J2112(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (-r1(y1, y2, y3) + (1 + 8 * ((
            - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + xlogy((-4)* ((-1) + nu)* ((
            - 1) + 2 * nu)* (x3 + y3), r2(y1, y2, y3) + x3 + y3)))
        end

        function J2113(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* ((x1 +
            (-1)* y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* (-((-1) +
            nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
            y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
            ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +
            xlogy(-x2 + y2, r1(y1, y2, y3) + x3 - y3)))
        end

        function J2123(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* (x1 + (
            -1)* y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (-((-1) +
            nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
            y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
            ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +
            xlogy(-x1 + y1, r1(y1, y2, y3) + x3 - y3)))
        end

        function J3112(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
            x3* (x2 - y2)* y3* (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
            -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)* atan((x1 - y1)
            * (x2 - y2)^ (-1)) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)*
            atan(r2(y1, y2, y3)* (-x1 + y1), (x2 - y2)* (x3 + y3)) + xlogy((-4)* ((-1) +
              nu)* ((-1) + 2 * nu)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +xlogy(x3 - y3, r1(y1, y2, y3) +
                x2 - y2) + xlogy(-x3 - 7 * y3 - 8 * nu^ 2 * (x3 + y3) + 8 * nu* (
                  x3 + 2 * y3), r2(y1, y2, y3) + x2 - y2)))
        end

        function J3113(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (r1(y1, y2, y3) + ((-1) - 8 * ((
            - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + 2 * ((-3) + 4 * nu)* x3*
            real(acoth(complex(r2(y1, y2, y3)^ (-1)* (x3 + y3)))) + xlogy(2 * (3 * x3 + 2 * y3 - 6 * nu* (x3 + y3) +
              4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3)))
        end

        function J3123(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
            -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)^ (-1)
            * (x2 - y2)) + 4 * ((-1) + 2 * nu)* (nu* x3 + ((-1) + nu)* y3)* atan(
              r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) + xlogy(x1 - y1, r1(y1, y2, y3) + x2 +
                (-1)* y2) + xlogy(-(1 + 8 * ((-1) + nu)* nu)* (x1 - y1), r2(y1, y2, y3) + x2 + (
                  -1)* y2)))
        end

        function J1212(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (-r1(y1, y2, y3) + (1 + 8 * ((
            - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + xlogy((-4)* ((-1) + nu)* ((
            - 1) + 2 * nu)* (x3 + y3), r2(y1, y2, y3) + x3 + y3)))
        end

        function J1213(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* ((x1 +
            (-1)* y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* (-((-1) +
            nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
            y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
            ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +
            xlogy(-x2 + y2, r1(y1, y2, y3) + x3 - y3)))
        end

        function J1223(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* (x1 + (
            -1)* y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (-((-1) +
            nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
            y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
            ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +
            xlogy(-x1 + y1, r1(y1, y2, y3) + x3 - y3)))
        end

        function J2212(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x1 - y1)* (x2 - y2)* y3* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (
            -1) -4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)* (
              x2 - y2)^ (-1)) -x3* atan(x3, x1 - y1) - 3 * x3*
            atan(3 * x3, x1 - y1) + 4 * nu* x3* atan(-nu* x3, x1 -
              y1) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan(r2(y1, y2, y3)* (-x2 + y2), (
                x1-y1)* (x3 + y3)) -4 * ((-1) + nu)* (x3 - y3)* atan(r1(y1, y2, y3)* (
                  x3 - y3), (x1 - y1)* (x2 - y2)) +3 * y3* atan((-3)* y3,
                    x1 - y1) - y3* atan(y3, x1 - y1) - 4 * nu* y3* atan(
                      nu* y3, x1 - y1) - 4 * ((-1) + nu)* (x3 + y3)* atan(r2(y1, y2, y3)* (x3 + y3), (
                        x1-y1)* (x2 - y2)) +xlogy((-4)* ((-1) + nu)* (x1 - y1),
                          r1(y1, y2, y3) + x2 - y2) + xlogy((-4)* ((-1) + nu)* (x1 - y1), r2(y1, y2, y3) + x2 -
                            y2) + xlogy(-((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x1 - y1) + xlogy(
                              (5 + 4 * nu* ((-3) + 2 * nu))* (x2 - y2), r2(y1, y2, y3) + x1 - y1)))
        end

        function J2213(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)* (
            x1 - y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* ((x2 + (-1)
            * y2)^ 2 + (x3 + y3)^ 2)^ (-1)* (x3* ((x3^ 2 + (x2 - y2)^ 2)* (
            x3^ 2 + (x1 - y1)^ 2 + (x2 - y2)^ 2) +x3* (3 * x3^ 2 + (x1 -
            y1)^ 2 + 2 * (x2 - y2)^ 2)* y3 + 3 * x3^ 2 * y3^ 2 + x3* y3^ 3) -(
            (-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)
            ^ 2) +((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)* y3* (2 * x3 + y3)* ((x2 - y2)
            ^ 2 + (x3 + y3)^ 2)) +2 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)* atan((
              x1-y1)^ (-1)* (x2 - y2)) +x2* atan(-x2, x1 - y1)
            -3 * x2* atan(3 * x2, x1 - y1) + 4 * nu* x2* atan(-nu* x2,
              x1 - y1) + 3 * y2* atan((-3)* y2, x1 - y1) - y2* atan(
                y2, x1 - y1) - 4 * nu* y2* atan(nu* y2, x1 - y1) + 2 * ((-1) +
            2 * nu)* (x2 - y2)* atan(r1(y1, y2, y3)* (-x2 + y2), (x1 - y1)* (x3 +
              (-1)* y3)) +2 * (1 - 2 * nu)^ 2 * (x2 - y2)* atan(r2(y1, y2, y3)* (-x2 +
                y2), (x1 - y1)* (x3 + y3)) +xlogy((-2)* x3, r2(y1, y2, y3) - x1 + y1) + xlogy((
            -1)* ((-3) + 4 * nu)* (x1 - y1), r1(y1, y2, y3) + x3 - y3) +xlogy(-(3 + (
                  -6)* nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +xlogy(-((-3) + 4 *
                    nu)* (x3 - y3), r1(y1, y2, y3) + x1 - y1) +xlogy(-(5 + 4 * nu* ((-3) + 2 *
                      nu))* (x3 + y3), r2(y1, y2, y3) + x1 - y1)))
        end

        function J2223(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* ((x1 +
            (-1)* y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* (-((-1) +
            nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
            y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
            ) +x1* atan(-x1, x2 - y2) - 3 * x1* atan(3 * x1, x2 -
              y2) + 4 * nu* x1* atan(-nu* x1, x2 - y2) - 4 * ((-1) + nu)* (
            x1 - y1)* atan(r1(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 - y3)
            ) + 4 * ((-1) + nu)* (x1 - y1)* atan(r2(y1, y2, y3)* (x1 - y1), (x2 -
              y2)* (x3 + y3)) +3 * y1* atan((-3)* y1, x2 - y2) - y1* atan(
                y1, x2 - y2) - 4 * nu* y1* atan(nu* y1, x2 - y2) + xlogy((-1)
            * ((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x3 - y3) + xlogy(-(3
                  -6 * nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +xlogy((-4)* ((-1) + nu)* (
                    x3 - y3), r1(y1, y2, y3) + x2 - y2) +xlogy(4 * ((-1) + nu)* (x3 + y3), r2(y1, y2, y3) + x2 + (
                      -1)* y2)))
        end

        function J3212(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
            x3* (x1 - y1)* y3* (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (
            -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)* atan((x1 - y1)
            ^ (-1)* (x2 - y2)) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)*
            atan(r2(y1, y2, y3)* (-x2 + y2), (x1 - y1)* (x3 + y3)) + xlogy((-4)* ((-1) +
              nu)* ((-1) + 2 * nu)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +xlogy(x3 - y3, r1(y1, y2, y3) +
                x1 - y1) + xlogy(-x3 - 7 * y3 - 8 * nu^ 2 * (x3 + y3) + 8 * nu* (
                  x3 + 2 * y3), r2(y1, y2, y3) + x1 - y1)))
        end

        function J3213(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x1 - y1)* (x2 - y2)* y3* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (
            -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)* (x2 + (
              -1)* y2)^ (-1)) +4 * ((-1) + 2 * nu)* (nu* x3 + ((-1) + nu)* y3)* atan(
                r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) + xlogy(x2 - y2, r1(y1, y2, y3) + x1 +
                  (-1)* y1) + xlogy(-(1 + 8 * ((-1) + nu)* nu)* (x2 - y2), r2(y1, y2, y3) + x1 + (
                    -1)* y1)))
        end

        function J3223(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (r1(y1, y2, y3) + ((-1) - 8 * ((
            - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + 2 * ((-3) + 4 * nu)* x3*
            real(acoth(complex(r2(y1, y2, y3)^ (-1)* (x3 + y3)))) + xlogy(2 * (3 * x3 + 2 * y3 - 6 * nu* (x3 + y3) +
              4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3)))
        end

        function J1312(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x2 - y2)* y3* (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (-1) + (
            -4)* ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)* atan((x1 - y1)* (
              x2 - y2)^ (-1)) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)*
            atan(r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) + xlogy(4 * ((-1) + nu)
            * ((-1) + 2 * nu)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) + xlogy(x3 - y3, r1(y1, y2, y3) + x2 + (
              -1)* y2) +xlogy((7 + 8 * ((-2) + nu)* nu)* x3 + y3 + 8 * ((-1) + nu)* nu* y3,
                r2(y1, y2, y3) + x2 - y2)))
        end

        function J1313(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (r1(y1, y2, y3) + r2(y1, y2, y3)^ (-1)* ((7 +
            8 * ((-2) + nu)* nu)* r2(y1, y2, y3)^ 2 + 2 * x3* y3) +2 * ((-3) + 4 * nu)* x3* real(acoth(
              complex(r2(y1, y2, y3)^ (-1)* (x3 + y3)))) + xlogy(2 * ((-3)* x3 - 2 * y3 + 6 * nu* (x3 + y3)
                -4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3)))
        end

        function J1323(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
            x3* (x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)
            ^ 2)^ (-1) - 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 -
              y1)^ (-1)* (x2 - y2)) -4 * ((-1) + nu)* ((-3)* x3 - y3 + 2 *
            nu* (x3 + y3))* atan(r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) +
            xlogy(x1 - y1, r1(y1, y2, y3) + x2 - y2) + xlogy((7 + 8 * ((-2) + nu)* nu)* (x1 +
              (-1)* y1), r2(y1, y2, y3) + x2 - y2)))
        end

        function J2312(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x1 - y1)* y3* (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (-1) + (
            -4)* ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)* atan((x1 - y1)^ (
              -1)* (x2 - y2)) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)*
            atan(r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) + xlogy(4 * ((-1) + nu)
            * ((-1) + 2 * nu)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) + xlogy(x3 - y3, r1(y1, y2, y3) + x1 + (
              -1)* y1) +xlogy((7 + 8 * ((-2) + nu)* nu)* x3 + y3 + 8 * ((-1) + nu)* nu* y3,
                r2(y1, y2, y3) + x1 - y1)))
        end

        function J2313(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
            x3* (x1 - y1)* (x2 - y2)* y3* ((x2 - y2)^ 2 + (x3 + y3)
            ^ 2)^ (-1) - 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 -
              y1)* (x2 - y2)^ (-1)) -4 * ((-1) + nu)* ((-3)* x3 - y3 + 2 *
            nu* (x3 + y3))* atan(r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) +
            xlogy(x2 - y2, r1(y1, y2, y3) + x1 - y1) + xlogy((7 + 8 * ((-2) + nu)* nu)* (x2 +
              (-1)* y2), r2(y1, y2, y3) + x1 - y1)))
        end

        function J2323(y1::R, y2::R, y3::R) where R
            ((-1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (r1(y1, y2, y3) + r2(y1, y2, y3)^ (-1)* ((7 +
            8 * ((-2) + nu)* nu)* r2(y1, y2, y3)^ 2 + 2 * x3* y3) +2 * ((-3) + 4 * nu)* x3* real(acoth(
              complex(r2(y1, y2, y3)^ (-1)* (x3 + y3)))) + xlogy(2 * ((-3)* x3 - 2 * y3 + 6 * nu* (x3 + y3)
                -4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3)))
        end

        function J3312(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
            -1)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (-1)* ((x1 - y1)^ 2 + (x2 + (
            -1)* y2)^ 2 + 2 * (x3 + y3)^ 2) -3 * x3* atan(3 * x3, x1 - y1)
            -5 * x3* atan(5 * x3, x2 - y2) + 12 * nu* x3* atan((-3)* nu* x3, x2 + (
              -1)* y2) +4 * nu* x3* atan(-nu* x3, x1 - y1) - 8 * nu^ 2 *
            x3* atan(nu^ 2 * x3, x2 - y2) + 3 * y3* atan((-3)* y3, x1 -
              y1) - 5 * y3* atan(5 * y3, x2 - y2) + 12 * nu* y3* atan((-3)*
                nu* y3, x2 - y2) - 4 * nu* y3* atan(nu* y3, x1 - y1) - 8 *
            nu^ 2 * y3* atan(nu^ 2 * y3, x2 - y2) + 2 * ((-1) + 2 * nu)* (x3 + (-1)
            * y3)* atan(r1(y1, y2, y3)* (-x3 + y3), (x1 - y1)* (x2 - y2)) + 2 * (
            1 - 2 * nu)^ 2 * (x3 + y3)* atan(r2(y1, y2, y3)* (x3 + y3), (x1 - y1)* (x2 + (-1)
            * y2)) +xlogy(-((-3) + 4 * nu)* (x1 - y1), r1(y1, y2, y3) + x2 - y2) +
            xlogy((5 + 4 * nu* ((-3) + 2 * nu))* (x1 - y1), r2(y1, y2, y3) + x2 - y2) +
            xlogy(-((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x1 - y1) + xlogy((5 +
              4 * nu* ((-3) + 2 * nu))* (x2 - y2), r2(y1, y2, y3) + x1 - y1)))
        end

        function J3313(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x1 - y1)* y3* (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (-1) + 5 *
            x2* atan((-5)* x2, x1 - y1) - 3 * x2* atan(3 * x2, x1 - y1)
            +4 * nu* x2* atan(-nu* x2, x1 - y1) - 12 * nu* x2* atan(
              3 * nu* x2, x1 - y1) + 8 * nu^ 2 * x2* atan(-nu^ 2 * x2, x1 + (-1)
            * y1) - 4 * ((-1) + nu)* (x2 - y2)* atan(r1(y1, y2, y3)* (x2 - y2), (x1 +
                (-1)* y1)* (x3 - y3)) -8 * ((-1) + nu)^ 2 * (x2 - y2)*
            atan(r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) + 3 * y2* atan((-3)
            * y2, x1 - y1) - 5 * y2* atan(5 * y2, x1 - y1) + 12 * nu* y2*
            atan((-3)* nu* y2, x1 - y1) - 4 * nu* y2* atan(nu* y2, x1 + (-1)
            * y1) - 8 * nu^ 2 * y2* atan(nu^ 2 * y2, x1 - y1) + xlogy((-4)*
              x3, r2(y1, y2, y3) - x1 + y1) + xlogy((-4)* ((-1) + nu)* (x1 - y1), r1(y1, y2, y3) + x3 + (-1)
            * y3) + xlogy((-8)* ((-1) + nu)^ 2 * (x1 - y1), r2(y1, y2, y3) + x3 + y3) + xlogy((-1)
            * ((-3) + 4 * nu)* (x3 - y3), r1(y1, y2, y3) + x1 - y1) + xlogy((-7)* x3
                -5 * y3 + 12 * nu* (x3 + y3) - 8 * nu^ 2 * (x3 + y3), r2(y1, y2, y3) + x1 - y1)))
        end

        function J3323(y1::R, y2::R, y3::R) where R
            ((1 / 16)* (1 - nu)^ (-1)* 1 / π * G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
            x2 - y2)* y3* (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (-1) + 5 *
            x1* atan((-5)* x1, x2 - y2) - 3 * x1* atan(3 * x1, x2 - y2)
            +4 * nu* x1* atan(-nu* x1, x2 - y2) - 12 * nu* x1* atan(
              3 * nu* x1, x2 - y2) + 8 * nu^ 2 * x1* atan(-nu^ 2 * x1, x2 + (-1)
            * y2) - 4 * ((-1) + nu)* (x1 - y1)* atan(r1(y1, y2, y3)* (x1 - y1), (x2 +
                (-1)* y2)* (x3 - y3)) -8 * ((-1) + nu)^ 2 * (x1 - y1)*
            atan(r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) + 3 * y1* atan((-3)
            * y1, x2 - y2) - 5 * y1* atan(5 * y1, x2 - y2) + 12 * nu* y1*
            atan((-3)* nu* y1, x2 - y2) - 4 * nu* y1* atan(nu* y1, x2 + (-1)
            * y2) - 8 * nu^ 2 * y1* atan(nu^ 2 * y1, x2 - y2) + xlogy((-4)*
              x3, r2(y1, y2, y3) - x2 + y2) + xlogy((-4)* ((-1) + nu)* (x2 - y2), r1(y1, y2, y3) + x3 + (-1)
            * y3) + xlogy((-8)* ((-1) + nu)^ 2 * (x2 - y2), r2(y1, y2, y3) + x3 + y3) + xlogy((-1)
            * ((-3) + 4 * nu)* (x3 - y3), r1(y1, y2, y3) + x2 - y2) + xlogy((-7)* x3
                -5 * y3 + 12 * nu* (x3 + y3) - 8 * nu^ 2 * (x3 + y3), r2(y1, y2, y3) + x2 - y2)))
        end

        function IU1(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J1123(y1, y2, y3)
            +2 * G * epsv12p * (J1223(y1, y2, y3) + J1113(y1, y2, y3))
            +2 * G * epsv13p * (J1323(y1, y2, y3) + J1112(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J1213(y1, y2, y3)
            +2 * G * epsv23p * (J1212(y1, y2, y3) + J1313(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J1312(y1, y2, y3))
        end

        function IU2(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J2123(y1, y2, y3)
            +2 * G * epsv12p * (J2223(y1, y2, y3) + J2113(y1, y2, y3))
            +2 * G * epsv13p * (J2323(y1, y2, y3) + J2112(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J2213(y1, y2, y3)
            +2 * G * epsv23p * (J2212(y1, y2, y3) + J2313(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J2312(y1, y2, y3))
        end

        function IU3(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J3123(y1, y2, y3)
            +2 * G * epsv12p * (J3223(y1, y2, y3) + J3113(y1, y2, y3))
            +2 * G * epsv13p * (J3323(y1, y2, y3) + J3112(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J3213(y1, y2, y3)
            +2 * G * epsv23p * (J3212(y1, y2, y3) + J3313(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J3312(y1, y2, y3))
        end

        u[1] =
            (IU1(L, T / 2, q3 + W) - IU1(L, -T / 2, q3 + W) + IU1(L, -T / 2, q3) - IU1(L, T / 2, q3)
            -IU1(zero(R), T / 2, q3 + W) + IU1(zero(R), -T / 2, q3 + W) - IU1(zero(R), -T / 2, q3) + IU1(zero(R), T / 2, q3))
        u[2] =
            (IU2(L, T / 2, q3 + W) - IU2(L, -T / 2, q3 + W) + IU2(L, -T / 2, q3) - IU2(L, T / 2, q3)
            -IU2(zero(R), T / 2, q3 + W) + IU2(zero(R), -T / 2, q3 + W) - IU2(zero(R), -T / 2, q3) + IU2(zero(R), T / 2, q3))
        u[3] =
            (IU3(L, T / 2, q3 + W) - IU3(L, -T / 2, q3 + W) + IU3(L, -T / 2, q3) - IU3(L, T / 2, q3)
            -IU3(zero(R), T / 2, q3 + W) + IU3(zero(R), -T / 2, q3 + W) - IU3(zero(R), -T / 2, q3) + IU3(zero(R), T / 2, q3))

        t1 = u[1] * cosd(theta) - u[2] * sind(theta)
        u[2] = u[1] * sind(theta) + u[2] * cosd(theta)
        u[1] = t1

        return nothing
    end
end

function sbarbot_strain_quad8!(ϵ::AbstractVector,
    x1::R, x2::R, x3::R, q1::R, q2::R, q3::R,
    L::R, T::R, W::R, theta::R,
    epsv11p::R, epsv12p::R, epsv13p::R, epsv22p::R, epsv23p::R, epsv33p::R,
    G::R, nu::R,
    ) where R

    lambda = G * 2 * nu / (1 - 2 * nu)
    epsvkk = epsv11p + epsv22p + epsv33p

    t1 = (x1 - q1) * cosd(theta) + (x2 - q2) * sind(theta)
    x2 = -(x1 - q1) * sind(theta) + (x2 - q2) * cosd(theta)
    x1 = t1

    # https://github com/JuliaLang/julia/issues/15276
    let lambda=lambda, epsvkk=epsvkk, x1=x1, x2=x2, x3=x3, q1=q1, q2=q2, q3=q3, L=L, T=T, W=W, theta=theta, epsv11p=epsv11p, epsv12p=epsv12p, epsv13p=epsv13p, epsv22p=epsv22p, epsv23p=epsv23p, epsv33p=epsv33p, G=G, nu=nu

        function r1(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 - y3)
        end

        function r2(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 + y3)
        end

        function J1112d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x3 ^ 2  * (x3 ^ 2 + (x1 +
            (-1) * y1) ^ 2) ^ (-1) + 9  * x3 ^ 2  * (9  * x3 ^ 2 + (x1 - y1) ^ 2) ^ (-1) +
            4  * nu ^ 2  * x3 ^ 2  * (nu ^ 2  * x3 ^ 2 + (x1 - y1) ^ 2) ^ (-1) - 4  * ((-1)
            +nu) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) - 4  * ((-1) + nu) * lr1 ^ (
            -1) * (x1 - y1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) - 4  * ((
                        - 1) +nu) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) - 4  * ((-1) + nu) *
            lr2 ^ (-1) * (x1 - y1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) + (
            -1) * ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * (lr1 + x2 - y2) ^ (
            -1) +(5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x2 + (
            -1) * y2) ^ (-1) + 4  * ((-1) + nu) * lr1 * (x2 - y2) * ((x1 - y1) ^ 2 +
            (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1)
                         * (x3 - y3) ^ 2 - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * (
            x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 + 4  * ((-1) + nu) * (
            (-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) * (x3 + y3) + y3 ^ 2  * ((x1 - y1) ^ 2 + y3 ^ 2) ^ (-1) + 9  * y3 ^ 2  * ((x1 +
            (-1) * y1) ^ 2 + 9  * y3 ^ 2) ^ (-1) + 4  * nu ^ 2  * y3 ^ 2  * ((x1 - y1) ^ 2 +
            nu ^ 2  * y3 ^ 2) ^ (-1) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) ^ 2  * (x2 + (-1)
                         * y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * x3 * (
            x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) +
            nu) * ((-1) + 2  * nu) * lr2 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) *
            (x2 - y2) * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + (
            -4) * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * ((
                x1- y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 * (x2 - y2)
                         * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1)
                         ^ 2  * (x2 - y2) * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1)
                         ^ 2  * (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 +
            (-1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + xlogy(3 - 4  *
                nu, lr1 + x2 - y2) + xlogy(5 + 4  * nu * ((-3) + 2  * nu), lr2 + x2 - y2)))
        end

        function J1112d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * ((-4) * ((-1) + nu) *
            lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) ^ 2 - 4  * ((-1) + nu)
                         * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) ^ 2 - ((-3) +
            4  * nu) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1) - ((-3) + 4  * nu) *
            lr1 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr1 + x2 - y2) ^ (-1) + (5 +
            4  * nu * ((-3) + 2  * nu)) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1) + (5 + 4  *
            nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x2 +
            (-1) * y2) ^ (-1) + 4  * ((-1) + nu) * lr1 * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) ^ 2 - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1)
                         * (x2 - y2) ^ 2  * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((
                x2- y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 - 4  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x3 + y3) + 2  * lr2 ^ (-1) * x3 * (x1 - y1) * y3 * ((x1 + (
            -1) * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 +
            (-1) * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) ^ 2  *
            ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu)
                         * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) ^ 2  * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-1) + 4  * ((-1) + nu) * lr2 * (x1 - y1) * (x3 + y3) ^ 2  * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * (
            (-1) + nu) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) ^ 2  * (x3 + y3) ^ 2  * (
            (x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-1) - 2  * x3 * (x1 - y1) * (x2 - y2) ^ 2  * y3 * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-3 / 2) + xlogy(4 - 4  * nu, lr1 + x1 - y1) + xlogy(4 - 4  * nu,
                lr2 + x1 - y1)))
        end

        function J1112d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x3 * (x3 ^ 2 + (x1 + (
            -1) * y1) ^ 2) ^ (-1) * (-x1 + y1) + 9  * x3 * (9  * x3 ^ 2 + (x1 - y1)
                         ^ 2) ^ (-1) * (-x1 + y1) + 4  * nu ^ 2  * x3 * (nu ^ 2  * x3 ^ 2 + (x1 -
            y1) ^ 2) ^ (-1) * (-x1 + y1) - 4  * ((-1) + nu) * lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) * (x3 - y3) - ((-3) + 4  * nu) *
            lr1 ^ (-1) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + (
            -4) * ((-1) + nu) * lr1 * (x1 - y1) * (x2 - y2) * ((x1 - y1)
                         ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) * (
            x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 3 - 4  * ((-1) + nu)
                         * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 + y3) + (5 + 4  *
            nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1)
                         * (x3 + y3) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (
            x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * (x1 + (-1)
                         * y1) * (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * (
            (-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            x3 + y3) ^ 3  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) *
            lr2 * (x1 - y1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu)
                         * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (x3 + y3) ^ 3  * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)
            -2  * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1) ^ (
                -1) * (x2 - y2)) -atan(x3, x1 - y1) - 3  * atan(3  * x3,
                    x1 - y1) + 4  * nu * atan(-nu * x3, x1 - y1) + 4  * ((-1) + nu) *
            ((-1) + 2  * nu) * atan(lr2 * (-x1 + y1), (x2 - y2) * (x3 + y3)) + (4 + (
            -4) * nu) * atan(lr1 * (x3 - y3), (x1 - y1) * (x2 - y2)) + (
            4 - 4  * nu) * atan(lr2 * (x3 + y3), (x1 - y1) * (x2 - y2))))
        end

        function J1113d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x2 ^ 2  * (x2 ^ 2 + (x1 +
            (-1) * y1) ^ 2) ^ (-1) + 9  * x2 ^ 2  * (9  * x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) +
            4  * nu ^ 2  * x2 ^ 2  * (nu ^ 2  * x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) + y2 ^ 2  * ((
                x1- y1) ^ 2 + y2 ^ 2) ^ (-1) + 9  * y2 ^ 2  * ((x1 - y1) ^ 2 + 9  * y2 ^ 2)
                         ^ (-1) + 4  * nu ^ 2  * y2 ^ 2  * ((x1 - y1) ^ 2 + nu ^ 2  * y2 ^ 2) ^ (-1)
            -4  * ((-1) + nu) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) - 4  * ((-1) + nu)
                         * lr1 ^ (-1) * (x1 - y1) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) +
            4  * ((-1) + nu) * lr1 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 -
            y3) -4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (
            x3 - y3) ^ 2) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) * lr1 ^ (-1)
                         * (x1 - y1) ^ 2  * (lr1 + x3 - y3) ^ (-1) + 4  * ((-1) + nu) * (lr2 + x1 + (
            -1) * y1) ^ (-1) * (x3 + y3) + 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * (
            lr2 + x1 - y1) ^ (-1) * (x3 + y3) - (3 - 6  * nu + 4  * nu ^ 2) * lr2 ^ (
            -1) * (x1 - y1) ^ 2  * (lr2 + x3 + y3) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) *
            y3 * (2  * x3 + y3) + 2  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * y3 * (2  * x3 + y3) + 2  * ((-1) + nu) * ((-1) + 2  * nu) *
            lr2 ^ (-2) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * y3 * (2  * x3 + y3) - 4  * ((-1) + nu) * lr2 * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((
                x2- y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (
            x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * ((-3) * nu * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (6  * nu * (x3 + y3) * (3  * (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) -4  * nu ^ 2  * (x3 + y3) * (3  * (x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2 + (x3 + y3) ^ 2) -2  * y3 * (3  * (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +xlogy(3 - 4  * nu, lr1 + x3 - y3) +
            xlogy((-3) + 6  * nu - 4  * nu ^ 2, lr2 + x3 + y3)))
        end

        function J1113d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x2 * (x2 ^ 2 + (x1 + (
            -1) * y1) ^ 2) ^ (-1) * (-x1 + y1) + 9  * x2 * (9  * x2 ^ 2 + (x1 - y1)
                         ^ 2) ^ (-1) * (-x1 + y1) + 4  * nu ^ 2  * x2 * (nu ^ 2  * x2 ^ 2 + (x1 -
            y1) ^ 2) ^ (-1) * (-x1 + y1) - 4  * ((-1) + nu) * lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) * (x3 - y3) - 4  * ((-1) + nu) * lr1 *
            (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) + (
            -4) * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 3  * ((x2 - y2) ^ 2 + (x3 -
            y3) ^ 2) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 + (
            -1) * y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) + 4  * ((-1) + nu) *
            lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 + y3) - (3 +
            (-6) * nu + 4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 +
            x3 + y3) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3)
            +2  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3)
            -4  * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 4  * ((
                        - 1) +nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 3  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) + 2  * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((
                x1- y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 - y3) * (
            2  * x3 + y3))) +4  * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-2) * (x2 - y2) * ((-3) * nu * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +atan(x1 - y1, (-1) * x2) - 3  *
            atan(3  * x2, x1 - y1) + 4  * nu * atan(-nu * x2, x1 - y1) + (
            4 - 4  * nu) * atan(lr1 * (x2 - y2), (x1 - y1) * (x3 - y3))
            +4  * ((-1) + nu) * atan(lr2 * (x2 - y2), (x1 - y1) * (x3 + y3))))
        end

        function J1113d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * lr1 *
            (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) - 4  * ((-1) +
            nu) * lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) ^ 2 - 4  * ((
                        - 1) +nu) * lr1 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) ^ 2 - ((-3) + 4  * nu) * (x1 - y1) * (lr1 +
            x3 - y3) ^ (-1) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 - y1) * (
            x3 - y3) * (lr1 + x3 - y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 + 4  * ((
                        - 1) +nu) * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x3 + y3) ^ 2 - (3 + (
            -6) * nu + 4  * nu ^ 2) * (x1 - y1) * (lr2 + x3 + y3) ^ (-1) - (3 - 6  *
            nu + 4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) +
            2  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 * (x3 + y3) * (2  * x3 + y3) - 2  *
            lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (3  * nu * ((-3) + 2  * nu) * x3 ^ 2 + nu * ((-3) + 2  * nu) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2) +2  * (2 - 9  * nu + 6  * nu ^ 2) * x3 * y3 + (3 - 9  * nu + 6  *
            nu ^ 2) * y3 ^ 2) -4  * ((-1) + nu) * lr2 * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) ^ 2  * ((
                x2- y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 *
            (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +xlogy(4 - 4  * nu, lr1 + x1 - y1) + xlogy(4  * (
                (-1) + nu), lr2 + x1 - y1)))
        end

        function J1123d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (2  * ((-1) + nu) * ((
                        - 1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) + x1 * (x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (-x2 +
            y2) +9  * x1 * (9  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (-x2 + y2) + 4  *
            nu ^ 2  * x1 * (nu ^ 2  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (-x2 + y2) +
            2  * lr2 ^ (-1) * x3 * (-x1 + y1) * (lr2 - x2 + y2) ^ (-1) - ((-3)
            +4  * nu) * lr1 ^ (-1) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1) * (x3 + (-1)
                         * y3) -2  * ((-1) + 2  * nu) * lr1 * (x1 - y1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 + (-1)
                         * y3) ^ 2) ^ (-1) * (x3 - y3) - 2  * ((-1) + 2  * nu) * lr1 ^ (-1) * (x1 + (
            -1) * y1) ^ 3  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) + (-1)
                         * ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr1 + x3 + (
            -1) * y3) ^ (-1) - (5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 -
            y1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) - (3 - 6  * nu + 4  * nu ^ 2) *
            lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + 4  * ((-1) +
            nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3) - 2  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * y3 * (2  * x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            (-2) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((
                        (-1) + 2  * nu) * lr2 * (x3 + y3) + 2  * ((-1) + nu) * y3 * (2  * x3 + y3)) +2  * lr2 ^ (-1)
                         * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (-x1 ^ 2  * x3 +
            2  * x1 * x3 * y1 - x3 * y1 ^ 2 + 3  * x1 ^ 2  * y3 + 2  * x2 ^ 2  * y3 + 8  * x3 ^ 2  *
            y3 - 6  * x1 * y1 * y3 + 3  * y1 ^ 2  * y3 - 4  * x2 * y2 * y3 + 2  * y2 ^ 2  * y3 +
            12  * x3 * y3 ^ 2 + 4  * y3 ^ 3 + 4  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + 2  * (x3 + y3) ^ 2) -2  * nu * (x3 + y3) * (4  * (x1 - y1)
                         ^ 2 + 3  * ((x2 - y2) ^ 2 + 2  * (x3 + y3) ^ 2))) -4  * lr2 ^ (-1) * (x1 + (-1)
                         * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            (x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) * (x1 ^ 4  * y3 - 4  * x1 ^ 3  * y1 *
            y3 - 3  * nu * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
            -2  * x1 * y1 * y3 * ((x2 - y2) ^ 2 + 2  * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3)
                        )) +y3 * (2  * x3 ^ 4 + 7  * x3 ^ 3  * y3 + (y1 ^ 2 + y3 ^ 2) * (y1 ^ 2 + y2 ^ 2 + y3 ^ 2) +
            x3 * y3 * (6  * y1 ^ 2 + 3  * y2 ^ 2 + 5  * y3 ^ 2) + x3 ^ 2  * (4  * y1 ^ 2 + 2  * y2 ^ 2 + 9  *
            y3 ^ 2) +x2 ^ 2  * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3)) - 2  * x2 * y2 * (y1 ^ 2 + (
            x3 + y3) * (2  * x3 + y3))) +x1 ^ 2  * y3 * ((x2 - y2) ^ 2 + 2  * (3  * y1 ^ 2 + (
            x3 + y3) * (2  * x3 + y3)))) -4  * lr2 ^ (-1) * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 - y2) * ((x1 - y1) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) * (x1 ^ 4  * y3 - 4  * x1 ^ 3  * y1 * y3 - 3  * nu * (x3 + y3)
                         * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) *
            ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) - 2  * x1 * y1 * y3 * (
            (x2 - y2) ^ 2 + 2  * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +y3 * (2  * x3 ^ 4 + 7  *
            x3 ^ 3  * y3 + (y1 ^ 2 + y3 ^ 2) * (y1 ^ 2 + y2 ^ 2 + y3 ^ 2) + x3 * y3 * (6  * y1 ^ 2 + 3  *
            y2 ^ 2 + 5  * y3 ^ 2) +x3 ^ 2  * (4  * y1 ^ 2 + 2  * y2 ^ 2 + 9  * y3 ^ 2) + x2 ^ 2  * (y1 ^ 2 +
            (x3 + y3) * (2  * x3 + y3)) -2  * x2 * y2 * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +
            x1 ^ 2  * y3 * ((x2 - y2) ^ 2 + 2  * (3  * y1 ^ 2 + (x3 + y3) * (2  * x3 + y3)))) +(
            -2) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (x1 ^ 4  * (y3 - 3  * nu * (
            x3 + y3) +2  * nu ^ 2  * (x3 + y3)) -4  * x1 ^ 3  * y1 * (y3 - 3  * nu * (x3 + y3) +
            2  * nu ^ 2  * (x3 + y3)) -3  * nu * (y1 ^ 2 + (x3 + y3) ^ 2) * (x3 ^ 3 + 3  * x3 ^ 2  *
            y3 + y3 * (y1 ^ 2 + y2 ^ 2 - lr2 * y3 + y3 ^ 2) + x3 * (y1 ^ 2 + y2 ^ 2 - 2  *
            lr2 * y3 + 3  * y3 ^ 2)) +2  * nu ^ 2  * (y1 ^ 2 + (x3 + y3) ^ 2) * (x3 ^ 3 + 3  * x3 ^ 2  *
            y3 + y3 * (y1 ^ 2 + y2 ^ 2 - lr2 * y3 + y3 ^ 2) + x3 * (y1 ^ 2 + y2 ^ 2 - 2  *
            lr2 * y3 + 3  * y3 ^ 2)) +y3 * (2  * x3 ^ 4 + 7  * x3 ^ 3  * y3 + (y1 ^ 2 + y3 ^ 2) * (
            y1 ^ 2 + y2 ^ 2 + y3 ^ 2) +x3 * y3 * (6  * y1 ^ 2 + 3  * y2 ^ 2 + 5  * y3 ^ 2) + x3 ^ 2  * (
            4  * y1 ^ 2 + 2  * y2 ^ 2 + 9  * y3 ^ 2) -lr2 * (2  * x3 + y3) * (y1 ^ 2 + (x3 + y3)
                         ^ 2)) +2  * x2 * y2 * (3  * nu * (x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2) - 2  * nu ^ 2  *
            (x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2) - y3 * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3)))
            +x2 ^ 2  * ((-3) * nu * (x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2) + 2  * nu ^ 2  * (x3 + y3) *
            (y1 ^ 2 + (x3 + y3) ^ 2) + y3 * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +2  * x1 * y1 * (((
                        - 1) + nu) * ((-1) + 2  * nu) * lr2 * y3 * (2  * x3 + y3) + x2 ^ 2  * (-y3 + 3  * nu *
            (x3 + y3) - 2  * nu ^ 2  * (x3 + y3)) +2  * x2 * y2 * (y3 - 3  * nu * (x3 + y3) + 2  *
            nu ^ 2  * (x3 + y3)) +3  * nu * (x3 + y3) * (2  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) + (
            -2) * nu ^ 2  * (x3 + y3) * (2  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) - y3 * (2  *
            y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) * (2  * x3 + y3))) +x1 ^ 2  * (-((-1) + nu) * ((
                        - 1) + 2  * nu) * lr2 * y3 * (2  * x3 + y3) + 2  * x2 * y2 * (-y3 + 3  * nu * (x3 + y3)
            -2  * nu ^ 2  * (x3 + y3)) +x2 ^ 2  * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 +
            y3)) -3  * nu * (x3 + y3) * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) + 2  * nu ^ 2  * (
            x3 + y3) * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) + y3 * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 +
            y3) * (2  * x3 + y3)))) +2  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1)
                         * (x2 - y2) ^ (-1)) + atan(x2 - y2, (-1) * x1) - 3  * atan(3  *
                x1, x2 - y2) + 4  * nu * atan(-nu * x1, x2 - y2) + ((-2) + 4  *
            nu) * atan(lr1 * (-x1 + y1), (x2 - y2) * (x3 - y3)) + 2  * (1 + (
            -2) * nu) ^ 2  * atan(lr2 * (-x1 + y1), (x2 - y2) * (x3 + y3))))
        end

        function J1123d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x1 ^ 2  * (x1 ^ 2 + (x2 +
            (-1) * y2) ^ 2) ^ (-1) + 9  * x1 ^ 2  * (9  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) +
            4  * nu ^ 2  * x1 ^ 2  * (nu ^ 2  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) - 2  * ((-1)
            +nu) * ((-1) + 2  * nu) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) + y1 ^ 2  * (y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 9  * y1 ^ 2  * (9  *
            y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 4  * nu ^ 2  * y1 ^ 2  * (nu ^ 2  * y1 ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) + 2  * x3 * (lr2 - x2 + y2) ^ (-1) + 2  * lr2 ^ (-1) * x3 *
            (-x2 + y2) * (lr2 - x2 + y2) ^ (-1) - ((-3) + 4  * nu) * (lr1 + x2 + (
            -1) * y2) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x2 + (
            -1) * y2) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + 2  * ((-1) + 2  * nu) *
            lr1 * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            (x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) - 2  * ((
                        - 1) +2  * nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x1 - y1) ^ 2 + (x3 - y3)
                         ^ 2) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x2 -
            y2) ^ 2  * (lr1 + x3 - y3) ^ (-1) - (5 + 4  * nu * ((-3) + 2  * nu)) * (lr2 +
            x2 - y2) ^ (-1) * (x3 + y3) - (5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (
            -1) * (x2 - y2) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) - (3
            -6  * nu + 4  * nu ^ 2) * lr2 ^ (-1) * (x2 - y2) ^ 2  * (lr2 + x3 + y3) ^ (-1) + 4  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) *
            (x2 - y2) ^ 2  * y3 * (2  * x3 + y3) - 2  * ((-1) + nu) * ((-1) + 2  * nu) *
            lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2)
                         ^ 2  * y3 * (2  * x3 + y3) + 2  * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (((-1) + 2  * nu)
                         * lr2 * (x1 - y1) ^ 2  * (x3 + y3) - ((-1) + nu) * y3 * (2  * x3 + y3) * (
            (x1 - y1) ^ 2 + (x3 + y3) ^ 2)) -4  * lr2 ^ (-1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) ^ (-2) * (x2 - y2) ^ 2  * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * (x1 ^ 4  * y3 - 4  * x1 ^ 3  * y1 * y3 - 3  * nu * (x3 + y3) * (
            (x1 - y1) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) * ((x1 +
            (-1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) -2  * x1 * y1 * y3 * ((x2 + (
            -1) * y2) ^ 2 + 2  * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +y3 * (2  * x3 ^ 4 + 7  *
            x3 ^ 3  * y3 + (y1 ^ 2 + y3 ^ 2) * (y1 ^ 2 + y2 ^ 2 + y3 ^ 2) + x3 * y3 * (6  * y1 ^ 2 + 3  *
            y2 ^ 2 + 5  * y3 ^ 2) +x3 ^ 2  * (4  * y1 ^ 2 + 2  * y2 ^ 2 + 9  * y3 ^ 2) + x2 ^ 2  * (y1 ^ 2 +
            (x3 + y3) * (2  * x3 + y3)) -2  * x2 * y2 * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +
            x1 ^ 2  * y3 * ((x2 - y2) ^ 2 + 2  * (3  * y1 ^ 2 + (x3 + y3) * (2  * x3 + y3)))) +(
            -2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  *
            ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (x1 ^ 4  * (y3 - 3  * nu * (x3 + y3) + 2  *
            nu ^ 2  * (x3 + y3)) -4  * x1 ^ 3  * y1 * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (
            x3 + y3)) -3  * nu * (y1 ^ 2 + (x3 + y3) ^ 2) * (x3 ^ 3 + 3  * x3 ^ 2  * y3 + y3 * (
            y1 ^ 2 + y2 ^ 2 - lr2 * y3 + y3 ^ 2) +x3 * (y1 ^ 2 + y2 ^ 2 - 2  * lr2 * y3 + 3  *
            y3 ^ 2)) +2  * nu ^ 2  * (y1 ^ 2 + (x3 + y3) ^ 2) * (x3 ^ 3 + 3  * x3 ^ 2  * y3 + y3 * (
            y1 ^ 2 + y2 ^ 2 - lr2 * y3 + y3 ^ 2) +x3 * (y1 ^ 2 + y2 ^ 2 - 2  * lr2 * y3 + 3  *
            y3 ^ 2)) +y3 * (2  * x3 ^ 4 + 7  * x3 ^ 3  * y3 + (y1 ^ 2 + y3 ^ 2) * (y1 ^ 2 + y2 ^ 2 +
            y3 ^ 2) +x3 * y3 * (6  * y1 ^ 2 + 3  * y2 ^ 2 + 5  * y3 ^ 2) + x3 ^ 2  * (4  * y1 ^ 2 + 2  *
            y2 ^ 2 + 9  * y3 ^ 2) -lr2 * (2  * x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2)) +2  * x2 *
            y2 * (3  * nu * (x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2) - 2  * nu ^ 2  * (x3 + y3) * (
            y1 ^ 2 + (x3 + y3) ^ 2) -y3 * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +x2 ^ 2  * ((
                        - 3) * nu * (x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2) + 2  * nu ^ 2  * (x3 + y3) * (y1 ^ 2 + (
            x3 + y3) ^ 2) +y3 * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +2  * x1 * y1 * (((-1) + nu)
                         * ((-1) + 2  * nu) * lr2 * y3 * (2  * x3 + y3) + x2 ^ 2  * (-y3 + 3  * nu * (x3 + y3)
            -2  * nu ^ 2  * (x3 + y3)) +2  * x2 * y2 * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  *
            (x3 + y3)) +3  * nu * (x3 + y3) * (2  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) - 2  *
            nu ^ 2  * (x3 + y3) * (2  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) - y3 * (2  *
            y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) * (2  * x3 + y3))) +x1 ^ 2  * (-((-1) + nu) * ((
                        - 1) + 2  * nu) * lr2 * y3 * (2  * x3 + y3) + 2  * x2 * y2 * (-y3 + 3  * nu * (x3 + y3)
            -2  * nu ^ 2  * (x3 + y3)) +x2 ^ 2  * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 +
            y3)) -3  * nu * (x3 + y3) * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) + 2  * nu ^ 2  * (
            x3 + y3) * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) + y3 * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 +
            y3) * (2  * x3 + y3)))) +2  * lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (-x3 *
            y1 ^ 2  * y2 ^ 2 + 2  * x3 ^ 4  * y3 + 4  * x3 ^ 2  * y1 ^ 2  * y3 + y1 ^ 4  * y3 + 6  * x3 ^ 2  *
            y2 ^ 2  * y3 + 2  * y1 ^ 2  * y2 ^ 2  * y3 + 7  * x3 ^ 3  * y3 ^ 2 + 6  * x3 * y1 ^ 2  * y3 ^ 2 +
            9  * x3 * y2 ^ 2  * y3 ^ 2 + 9  * x3 ^ 2  * y3 ^ 3 + 2  * y1 ^ 2  * y3 ^ 3 + 3  * y2 ^ 2  *
            y3 ^ 3 + 5  * x3 * y3 ^ 4 + y3 ^ 5 - nu * (x3 + y3) * (3  * x3 ^ 4 + 6  * x3 ^ 2  *
            y1 ^ 2 + 3  * y1 ^ 4 + 9  * x3 ^ 2  * y2 ^ 2 + 5  * y1 ^ 2  * y2 ^ 2 + 6  * x3 * (2  * (x3 ^ 2 +
            y1 ^ 2) +3  * y2 ^ 2) * y3 + 3  * (6  * x3 ^ 2 + 2  * y1 ^ 2 + 3  * y2 ^ 2) * y3 ^ 2 + 12  *
            x3 * y3 ^ 3 + 3  * y3 ^ 4) +x1 ^ 4  * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)
                        ) -4  * x1 ^ 3  * y1 * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 2  *
            nu ^ 2  * (x3 + y3) * (x3 ^ 4 + y1 ^ 4 + 4  * x3 ^ 3  * y3 + 3  * y2 ^ 2  * y3 ^ 2 + y3 ^ 4 +
            y1 ^ 2  * (y2 ^ 2 + 2  * y3 ^ 2) + x3 ^ 2  * (2  * y1 ^ 2 + 3  * y2 ^ 2 + 6  * y3 ^ 2) + x3 * (
            4  * y1 ^ 2  * y3 + 6  * y2 ^ 2  * y3 + 4  * y3 ^ 3)) +x2 ^ 2  * (-x3 * y1 ^ 2 + 6  *
            x3 ^ 2  * y3 + 2  * y1 ^ 2  * y3 + 9  * x3 * y3 ^ 2 + 3  * y3 ^ 3 + 2  * nu ^ 2  * (x3 + y3) * (
            y1 ^ 2 + 3  * (x3 + y3) ^ 2) -nu * (x3 + y3) * (5  * y1 ^ 2 + 9  * (x3 + y3) ^ 2)) +
            2  * x2 * y2 * (x3 * y1 ^ 2 - 6  * x3 ^ 2  * y3 - 2  * y1 ^ 2  * y3 - 9  * x3 *
            y3 ^ 2 - 3  * y3 ^ 3 - 2  * nu ^ 2  * (x3 + y3) * (y1 ^ 2 + 3  * (x3 + y3) ^ 2) + nu *
            (x3 + y3) * (5  * y1 ^ 2 + 9  * (x3 + y3) ^ 2)) +2  * x1 * y1 * (x2 ^ 2  * x3 - 2  *
            x2 * x3 * y2 + x3 * y2 ^ 2 - 2  * x2 ^ 2  * y3 - 4  * x3 ^ 2  * y3 - 2  * y1 ^ 2  *
            y3 + 4  * x2 * y2 * y3 - 2  * y2 ^ 2  * y3 - 6  * x3 * y3 ^ 2 - 2  * y3 ^ 3
            -2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + 2  * (y1 ^ 2 + (x3 + y3) ^ 2)) + nu * (
            x3 + y3) * (5  * (x2 - y2) ^ 2 + 6  * (y1 ^ 2 + (x3 + y3) ^ 2))) +x1 ^ 2  * ((-1)
                         * x2 ^ 2  * x3 + 2  * x2 * x3 * y2 - x3 * y2 ^ 2 + 2  * x2 ^ 2  * y3 + 4  * x3 ^ 2  *
            y3 + 6  * y1 ^ 2  * y3 - 4  * x2 * y2 * y3 + 2  * y2 ^ 2  * y3 + 6  * x3 * y3 ^ 2 + 2  *
            y3 ^ 3 + 2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + 2  * (3  * y1 ^ 2 + (x3 + y3)
                         ^ 2)) -nu * (x3 + y3) * (5  * (x2 - y2) ^ 2 + 6  * (3  * y1 ^ 2 + (x3 + y3)
                         ^ 2)))) +xlogy(3 - 4  * nu, lr1 + x3 - y3) + xlogy((-3) + 6  * nu - 4  *
                nu ^ 2, lr2 + x3 + y3)))
        end

        function J1123d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (2  * ((-1) + 2  * nu) *
            lr1 * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) - ((
                        - 3) +4  * nu) * lr1 ^ (-1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) ^ 2 + (
            -2) * ((-1) + 2  * nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 + (-1)
                         * y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 - ((-3) + 4  * nu) * (x2 -
            y2) * (lr1 + x3 - y3) ^ (-1) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x2 + (-1)
                         * y2) * (x3 - y3) * (lr1 + x3 - y3) ^ (-1) - 2  * lr2 ^ (-1) * x3 * (
            lr2 - x2 + y2) ^ (-1) * (x3 + y3) - (5 + 4  * nu * ((-3) + 2  * nu)) *
            lr2 ^ (-1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) ^ 2 - (3 - 6  * nu + 4  *
            nu ^ 2) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) - (3 - 6  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) - 2  * (
            (-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * y3 * (x3 + y3) * (2  * x3 + y3) + 4  * ((-1) + nu) * (
            (-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) * y3 * (x3 + y3) * (2  * x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            2  * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (((-1) + 2  * nu) * lr2 *
            (x1 - y1) ^ 2 - 2  * ((-1) + nu) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) * (
            3  * x3 + 2  * y3))) +2  * lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (nu * ((
                        - 3) + 2  * nu) * x1 ^ 4 + 4  * (3 - 2  * nu) * nu * x1 ^ 3  * y1 - x3 ^ 2  *
            y1 ^ 2 + 4  * x2 ^ 2  * x3 * y3 + 8  * x3 ^ 3  * y3 + 6  * x3 * y1 ^ 2  * y3 - 8  * x2 *
            x3 * y2 * y3 + 4  * x3 * y2 ^ 2  * y3 + 3  * x2 ^ 2  * y3 ^ 2 + 21  * x3 ^ 2  * y3 ^ 2 + 5  *
            y1 ^ 2  * y3 ^ 2 - 6  * x2 * y2 * y3 ^ 2 + 3  * y2 ^ 2  * y3 ^ 2 + 18  * x3 * y3 ^ 3 + 5  *
            y3 ^ 4 + 2  * x1 * y1 * ((1 + 2  * (7 - 4  * nu) * nu) * x3 ^ 2 + (3 - 2  * nu) *
            nu * (2  * y1 ^ 2 + (x2 - y2) ^ 2) - 2  * (3 + 2  * nu * ((-7) + 4  * nu)) *
            x3 * y3 + ((-5) + 2  * (7 - 4  * nu) * nu) * y3 ^ 2) +x1 ^ 2  * (((-1) + 2  * nu * ((
                        - 7) + 4  * nu)) * x3 ^ 2 + nu * ((-3) + 2  * nu) * (6  * y1 ^ 2 + (x2 - y2) ^ 2) +
            2  * (3 + 2  * nu * ((-7) + 4  * nu)) * x3 * y3 + (5 + 2  * nu * ((-7) + 4  * nu)) * y3 ^ 2)
            -nu * (y1 ^ 2 + 3  * (x3 + y3) ^ 2) * (3  * y1 ^ 2 + 3  * (x2 - y2) ^ 2 +
            5  * (x3 + y3) ^ 2) +2  * nu ^ 2  * (5  * x3 ^ 4 + y1 ^ 4 + y1 ^ 2  * y2 ^ 2 + 20  * x3 ^ 3  *
            y3 + 4  * y1 ^ 2  * y3 ^ 2 + 3  * y2 ^ 2  * y3 ^ 2 + 5  * y3 ^ 4 + x3 ^ 2  * (4  * y1 ^ 2 + 3  *
            y2 ^ 2 + 30  * y3 ^ 2) +x3 * (8  * y1 ^ 2  * y3 + 6  * y2 ^ 2  * y3 + 20  * y3 ^ 3) + x2 ^ 2  *
            (y1 ^ 2 + 3  * (x3 + y3) ^ 2) - 2  * x2 * y2 * (y1 ^ 2 + 3  * (x3 + y3) ^ 2))) -4  *
            lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2)
                         * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) * (x1 ^ 4  * y3 - 4  *
            x1 ^ 3  * y1 * y3 - 3  * nu * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) * ((
                x1- y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((
                    x1- y1) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) -2  * x1 * y1 * y3 * ((x2 - y2) ^ 2 + 2  * (y1 ^ 2 + (x3 + y3) *
            (2  * x3 + y3))) +y3 * (2  * x3 ^ 4 + 7  * x3 ^ 3  * y3 + (y1 ^ 2 + y3 ^ 2) * (y1 ^ 2 +
            y2 ^ 2 + y3 ^ 2) +x3 * y3 * (6  * y1 ^ 2 + 3  * y2 ^ 2 + 5  * y3 ^ 2) + x3 ^ 2  * (4  *
            y1 ^ 2 + 2  * y2 ^ 2 + 9  * y3 ^ 2) +x2 ^ 2  * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3)) - 2  *
            x2 * y2 * (y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +x1 ^ 2  * y3 * ((x2 - y2) ^ 2 +
            2  * (3  * y1 ^ 2 + (x3 + y3) * (2  * x3 + y3)))) -2  * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -3 / 2) * (x1 ^ 4  * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) - 4  *
            x1 ^ 3  * y1 * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) - 3  * nu * (
            y1 ^ 2 + (x3 + y3) ^ 2) * (x3 ^ 3 + 3  * x3 ^ 2  * y3 + y3 * (y1 ^ 2 + y2 ^ 2 - lr2 *
            y3 + y3 ^ 2) +x3 * (y1 ^ 2 + y2 ^ 2 - 2  * lr2 * y3 + 3  * y3 ^ 2)) +2  * nu ^ 2  * (
            y1 ^ 2 + (x3 + y3) ^ 2) * (x3 ^ 3 + 3  * x3 ^ 2  * y3 + y3 * (y1 ^ 2 + y2 ^ 2 - lr2 *
            y3 + y3 ^ 2) +x3 * (y1 ^ 2 + y2 ^ 2 - 2  * lr2 * y3 + 3  * y3 ^ 2)) +y3 * (2  * x3 ^ 4 +
            7  * x3 ^ 3  * y3 + (y1 ^ 2 + y3 ^ 2) * (y1 ^ 2 + y2 ^ 2 + y3 ^ 2) + x3 * y3 * (6  * y1 ^ 2 +
            3  * y2 ^ 2 + 5  * y3 ^ 2) +x3 ^ 2  * (4  * y1 ^ 2 + 2  * y2 ^ 2 + 9  * y3 ^ 2) - lr2 * (
            2  * x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2)) +2  * x2 * y2 * (3  * nu * (x3 + y3) * (y1 ^ 2 +
            (x3 + y3) ^ 2) -2  * nu ^ 2  * (x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2) - y3 * (
            y1 ^ 2 + (x3 + y3) * (2  * x3 + y3))) +x2 ^ 2  * ((-3) * nu * (x3 + y3) * (y1 ^ 2 + (x3 +
            y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * (y1 ^ 2 + (x3 + y3) ^ 2) + y3 * (y1 ^ 2 + (x3 + y3) *
            (2  * x3 + y3))) +2  * x1 * y1 * (((-1) + nu) * ((-1) + 2  * nu) * lr2 * y3 * (2  * x3 +
            y3) +x2 ^ 2  * (-y3 + 3  * nu * (x3 + y3) - 2  * nu ^ 2  * (x3 + y3)) + 2  * x2 *
            y2 * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 3  * nu * (x3 + y3) * (2  *
            y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) ^ 2) -2  * nu ^ 2  * (x3 + y3) * (2  * y1 ^ 2 + y2 ^ 2 +
            2  * (x3 + y3) ^ 2) -y3 * (2  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) * (2  * x3 + y3))) +
            x1 ^ 2  * (-((-1) + nu) * ((-1) + 2  * nu) * lr2 * y3 * (2  * x3 + y3) + 2  * x2 *
            y2 * (-y3 + 3  * nu * (x3 + y3) - 2  * nu ^ 2  * (x3 + y3)) + x2 ^ 2  * (y3
            -3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) -3  * nu * (x3 + y3) * (6  * y1 ^ 2 +
            y2 ^ 2 + 2  * (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3)
                         ^ 2) +y3 * (6  * y1 ^ 2 + y2 ^ 2 + 2  * (x3 + y3) * (2  * x3 + y3)))) +xlogy((-2), lr2 + (
                -1) * x2 + y2) +xlogy(3 - 4  * nu, lr1 + x2 - y2) + xlogy((-5) + 4  * (3
                    -2  * nu) * nu, lr2 + x2 - y2)))
        end

        function J2112d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((1 + 8  * ((-1) + nu) * nu) * lr2 ^ (
            -1) * (x1 - y1) + lr1 ^ (-1) * (-x1 + y1) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * x3 *
            (x1 - y1) * y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2)) * G ^ (-1))
        end

        function J2112d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((1 + 8  * ((-1) + nu) * nu) * lr2 ^ (
            -1) * (x2 - y2) + lr1 ^ (-1) * (-x2 + y2) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * x3 *
            (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2)) * G ^ (-1))
        end

        function J2112d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (-
            x3 + y3) -4  * ((-1) + nu) * ((-1) + 2  * nu) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + (
            -4) * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x3 + y3) ^ 2  * (lr2 + x3 + y3) ^ (
            -1) +lr2 ^ (-1) * (x3 - y3 - 8  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) +
            2  * x3 * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2) + xlogy((-4) * ((-1) + nu) * ((-1) + 2  * nu), lr2 + x3 + y3)))
        end

        function J2113d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * (-lr1 ^ (-1) * (x1 -
            y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu +
            4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (
            -1) -4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) - 4  * lr2 ^ (
            -1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 2  * (x1 + (-1)
                         * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            (x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu *
            (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  *
            x3) +(x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) *
            y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 +
            (-1) * (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 + (-1)
                         * y2) * ((-3) * nu * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3)
                        ))) * G ^ (-1))
        end

        function J2113d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (-lr1 ^ (-1) * (
            x2 - y2) ^ 2  * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x2 - y2) ^ 2  * (lr2 + x3 + y3) ^ (-1) + 2  * ((-1) + nu)
                         * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 * (
            2  * x3 + y3) -4  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-2) * (x2 - y2) ^ 2  * y3 * (2  * x3 + y3) + 2  * ((-1) + nu)
                         * ((-1) + 2  * nu) * lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) ^ 2  * y3 * (2  * x3 + y3) + 2  * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 *
            (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-2) * (x2 - y2) ^ 2  * ((-3) * nu * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (6  * nu * (x3 + y3) * ((x1 - y1) ^ 2 + 3  * (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) -4  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + 3  * (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) -2  * y3 * ((x1 - y1) ^ 2 + 3  * (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +log(lr2 + x3 + y3) + xlogy((-1), lr1 + x3 + (
                -1) * y3) +xlogy(2  * (1 - 2  * nu) * nu, lr2 + x3 + y3)))
        end

        function J2113d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((-x2 + y2) * (lr1 + x3 + (-1)
                         * y3) ^ (-1) - lr1 ^ (-1) * (x2 - y2) * (x3 - y3) * (lr1 + x3 + (
            -1) * y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) * y3 - ((-1) - 2  * nu + 4  *
            nu ^ 2) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * y3 * (x3 + y3) * (2  * x3 + y3) - 2  * lr2 ^ (-1)
                         * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (3  *
            nu * ((-3) + 2  * nu) * x3 ^ 2 + nu * ((-3) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) +2  * (2 - 9  * nu + 6  * nu ^ 2) * x3 * y3 + (3 - 9  * nu + 6  *
            nu ^ 2) * y3 ^ 2) +2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 +
            (-1) * (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 +
            y3)))) * G ^ (-1))
        end

        function J2123d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (-lr1 ^ (-1) * (
            x1 - y1) ^ 2  * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x3 + y3) ^ (-1) - 4  * ((-1) +
            nu) * ((-1) + 2  * nu) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-2) * y3 * (2  * x3 + y3) + 2  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 * (2  * x3 + y3) + 2  * ((-1) + nu) *
            ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * y3 * (2  * x3 + y3) + 2  * (x1 - y1) ^ 2  * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (
            x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * ((-3) * nu * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (6  * nu * (x3 + y3) * (3  * (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) -4  * nu ^ 2  * (x3 + y3) * (3  * (x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2 + (x3 + y3) ^ 2) -2  * y3 * (3  * (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +log(lr2 + x3 + y3) + xlogy((-1), lr1 + x3 -
                y3) + xlogy(2  * (1 - 2  * nu) * nu, lr2 + x3 + y3)))
        end

        function J2123d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * (-lr1 ^ (-1) * (x1 -
            y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu +
            4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (
            -1) -4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) - 4  * lr2 ^ (
            -1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 2  * (x1 + (-1)
                         * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            (x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu *
            (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  *
            x3) +(x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) *
            y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 +
            (-1) * (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 + (-1)
                         * y2) * ((-3) * nu * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3)
                        ))) * G ^ (-1))
        end

        function J2123d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((-x1 + y1) * (lr1 + x3 + (-1)
                         * y3) ^ (-1) - lr1 ^ (-1) * (x1 - y1) * (x3 - y3) * (lr1 + x3 + (
            -1) * y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 - ((-1) - 2  * nu + 4  *
            nu ^ 2) * (x1 - y1) * (lr2 + x3 + y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * y3 * (x3 + y3) * (2  * x3 + y3) - 2  * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (3  * nu * (
            (-3) + 2  * nu) * x3 ^ 2 + nu * ((-3) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) +2  * (2 - 9  * nu + 6  * nu ^ 2) * x3 * y3 + (3 - 9  * nu + 6  * nu ^ 2) *
            y3 ^ 2) +2  * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2)
                         * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * (
            (-2) * lr2 + 3  * x3) +(x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 + (
            -3) * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 + y3)))) * G ^ (
            -1))
        end

        function J3112d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * ((
                        - 1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) + lr1 ^ (-1) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1)
                         * (x3 - y3) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 + (-1)
                         * y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + 4  * lr2 ^ (-1) * x3 * (x1 -
            y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (
            -2) -4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 + (-1)
                         * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1)
                         * (x1 - y1) ^ 3  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * x3 *
            (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -3 / 2) -lr2 ^ (-1) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1) * (x3 +
            7  * y3 + 8  * nu ^ 2  * (x3 + y3) - 8  * nu * (x3 + 2  * y3)) +4  * ((-1) + nu) * ((-1) +
            2  * nu) * atan((x1 - y1) * (x2 - y2) ^ (-1)) + 4  * ((-1) + nu) * ((
                        - 1) +2  * nu) * atan(lr2 * (-x1 + y1), (x2 - y2) * (x3 + y3))))
        end

        function J3112d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * ((-4) * ((-1) + nu) * (
            (-1) + 2  * nu) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         ^ (-1) + (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + lr1 ^ (-1) * (x2 -
            y2) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 ^ (-1) * (x2 - y2) ^ 2  * (lr2 + x3 + y3) ^ (-1) + 4  * ((-1) + nu)
                         * ((-1) + 2  * nu) * lr2 * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  *
            ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((x1 + (
            -1) * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * lr2 ^ (-1) * x3 * y3 * (x3 + y3) * ((
                x1- y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * x3 * (x2 - y2) ^ 2  * y3 * (
            x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - lr2 ^ (-1) * (x2 - y2)
                         * (lr2 + x2 - y2) ^ (-1) * (x3 + 7  * y3 + 8  * nu ^ 2  * (x3 + y3) - 8  * nu * (
            x3 + 2  * y3)) +(lr2 + x2 - y2) ^ (-1) * (-x3 - 7  * y3 - 8  *
            nu ^ 2  * (x3 + y3) + 8  * nu * (x3 + 2  * y3)) +xlogy((-4) * ((-1) + nu) * ((-1) + 2  *
                nu), lr2 + x3 + y3)))
        end

        function J3112d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x2 + (
            -1) * y2) ^ (-1) * (x3 - y3) ^ 2 - 4  * ((-1) + nu) * ((-1) + 2  * nu) * (
            x2 - y2) * (lr2 + x3 + y3) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) *
            lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 4  * lr2 ^ (-1) *
            x3 * (x2 - y2) * y3 * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-2) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 - y1) ^ 2  * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (
            x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  *
            lr2 ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 2  * x3 * (x2 - y2) * y3 * (x3 + y3) ^ 2  * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) - lr2 ^ (-1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) * (x3 +
            7  * y3 + 8  * nu ^ 2  * (x3 + y3) - 8  * nu * (x3 + 2  * y3)) +log(lr1 + x2 - y2) +
            xlogy((-1) - 8  * ((-1) + nu) * nu, lr2 + x2 - y2)))
        end

        function J3113d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x1 - y1) + (
            -1) * (1 + 8  * ((-1) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) + 2  * lr2 ^ (-1) * (
            x1 - y1) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) +2  * x3 * (x1 - y1) * y3 * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 2  * ((-3) + 4  * nu) * x3 * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2)) * G ^ (-1))
        end

        function J3113d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x2 - y2) + (
            -1) * (1 + 8  * ((-1) + nu) * nu) * lr2 ^ (-1) * (x2 - y2) + 2  * lr2 ^ (-1) * (
            x2 - y2) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) +2  * x3 * (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 2  * ((-3) + 4  * nu) * x3 * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2)) * G ^ (-1))
        end

        function J3113d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x3 + (-1)
                         * y3) +lr2 ^ (-1) * (-x3 - 3  * y3 + 8  * nu * (x3 + y3) - 8  * nu ^ 2  * (
            x3 + y3)) +2  * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) +2  * lr2 ^ (-1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  *
            y3 - 6  * nu * (x3 + y3) + 4  * nu ^ 2  * (x3 + y3)) +2  * x3 * y3 * (x3 + y3) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 2  * ((-3) + 4  *
            nu) * x3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) ^ 2  * ((
                x1- y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + 2  * ((-3) + 4  *
            nu) * lr2 ^ (-1) * x3 * (1 - (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)) ^ (-1) + 2  * ((-3) + 4  * nu) * real(acoth(complex(lr2 ^ (-1)
                         * (x3 + y3)))) + xlogy(6 + 4  * nu * ((-3) + 2  * nu), lr2 + x3 + y3)))
        end

        function J3123d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x1 + (-1)
                         * y1) ^ 2  * (lr1 + x2 - y2) ^ (-1) - (1 + 8  * ((-1) + nu) * nu) * lr2 ^ (
            -1) * (x1 - y1) ^ 2  * (lr2 + x2 - y2) ^ (-1) - 4  * ((-1) + nu) * ((
                        - 1) +2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) * (x3 + y3) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) ^ 2  * (x2 - y2) *
            y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * x3 * (x2 + (-1)
                         * y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + 2  * nu) *
            lr2 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 +
            y3) * (nu * x3 + ((-1) + nu) * y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            4  * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * (nu * x3 + ((-1) + nu)
                         * y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1)
                         ^ 2  * (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 +
            (-1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + log(lr1 + x2 -
                y2) + xlogy((-1) - 8  * ((-1) + nu) * nu, lr2 + x2 - y2)))
        end

        function J3123d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((x1 - y1) * (lr1 + x2 + (-1)
                         * y2) ^ (-1) + lr1 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr1 + x2 -
            y2) ^ (-1) - (1 + 8  * ((-1) + nu) * nu) * (x1 - y1) * (lr2 + x2 -
            y2) ^ (-1) - (1 + 8  * ((-1) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) * (
            x2 - y2) * (lr2 + x2 - y2) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) +
            2  * lr2 ^ (-1) * x3 * (x1 - y1) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-1) - 4  * ((-1) + 2  * nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * (nu * x3 + ((-1) + nu) * y3) * ((x1 + (-1)
                         * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  *
            (x3 + y3) * (nu * x3 + ((-1) + nu) * y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) -2  * x3 * (x1 - y1) * (x2 - y2) ^ 2  * y3 * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2)) * G ^ (-1))
        end

        function J3123d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x1 + (-1)
                         * y1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) - (1 + 8  * ((-1) + nu)
                         * nu) * lr2 ^ (-1) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) + (
            -4) * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((
                x1- y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * (x1 - y1) * (x2 +
            (-1) * y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + 2  *
            nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x2 - y2) * (nu * x3 + ((-1) + nu) * y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 4  * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * (nu *
            x3 + ((-1) + nu) * y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (
            x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -3 / 2) +4  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1) ^ (-1) * (x2 + (
                -1) * y2)) +4  * nu * ((-1) + 2  * nu) * atan(lr2 * (x1 - y1), (x2 -
                    y2) * (x3 + y3))))
        end

        function J1212d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((1 + 8  * ((-1) + nu) * nu) * lr2 ^ (
            -1) * (x1 - y1) + lr1 ^ (-1) * (-x1 + y1) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * x3 *
            (x1 - y1) * y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2)) * G ^ (-1))
        end

        function J1212d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((1 + 8  * ((-1) + nu) * nu) * lr2 ^ (
            -1) * (x2 - y2) + lr1 ^ (-1) * (-x2 + y2) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * x3 *
            (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2)) * G ^ (-1))
        end

        function J1212d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (-
            x3 + y3) -4  * ((-1) + nu) * ((-1) + 2  * nu) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + (
            -4) * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x3 + y3) ^ 2  * (lr2 + x3 + y3) ^ (
            -1) +lr2 ^ (-1) * (x3 - y3 - 8  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) +
            2  * x3 * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2) + xlogy((-4) * ((-1) + nu) * ((-1) + 2  * nu), lr2 + x3 + y3)))
        end

        function J1213d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * (-lr1 ^ (-1) * (x1 -
            y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu +
            4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (
            -1) -4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) - 4  * lr2 ^ (
            -1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 2  * (x1 + (-1)
                         * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            (x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu *
            (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  *
            x3) +(x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) *
            y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 +
            (-1) * (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 + (-1)
                         * y2) * ((-3) * nu * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3)
                        ))) * G ^ (-1))
        end

        function J1213d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (-lr1 ^ (-1) * (
            x2 - y2) ^ 2  * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x2 - y2) ^ 2  * (lr2 + x3 + y3) ^ (-1) + 2  * ((-1) + nu)
                         * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 * (
            2  * x3 + y3) -4  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-2) * (x2 - y2) ^ 2  * y3 * (2  * x3 + y3) + 2  * ((-1) + nu)
                         * ((-1) + 2  * nu) * lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) ^ 2  * y3 * (2  * x3 + y3) + 2  * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 *
            (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-2) * (x2 - y2) ^ 2  * ((-3) * nu * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (6  * nu * (x3 + y3) * ((x1 - y1) ^ 2 + 3  * (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) -4  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + 3  * (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) -2  * y3 * ((x1 - y1) ^ 2 + 3  * (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +log(lr2 + x3 + y3) + xlogy((-1), lr1 + x3 + (
                -1) * y3) +xlogy(2  * (1 - 2  * nu) * nu, lr2 + x3 + y3)))
        end

        function J1213d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((-x2 + y2) * (lr1 + x3 + (-1)
                         * y3) ^ (-1) - lr1 ^ (-1) * (x2 - y2) * (x3 - y3) * (lr1 + x3 + (
            -1) * y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) * y3 - ((-1) - 2  * nu + 4  *
            nu ^ 2) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * y3 * (x3 + y3) * (2  * x3 + y3) - 2  * lr2 ^ (-1)
                         * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (3  *
            nu * ((-3) + 2  * nu) * x3 ^ 2 + nu * ((-3) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) +2  * (2 - 9  * nu + 6  * nu ^ 2) * x3 * y3 + (3 - 9  * nu + 6  *
            nu ^ 2) * y3 ^ 2) +2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 +
            (-1) * (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 +
            y3)))) * G ^ (-1))
        end

        function J1223d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (-lr1 ^ (-1) * (
            x1 - y1) ^ 2  * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x3 + y3) ^ (-1) - 4  * ((-1) +
            nu) * ((-1) + 2  * nu) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-2) * y3 * (2  * x3 + y3) + 2  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 * (2  * x3 + y3) + 2  * ((-1) + nu) *
            ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * y3 * (2  * x3 + y3) + 2  * (x1 - y1) ^ 2  * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (
            x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * ((-3) * nu * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (6  * nu * (x3 + y3) * (3  * (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) -4  * nu ^ 2  * (x3 + y3) * (3  * (x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2 + (x3 + y3) ^ 2) -2  * y3 * (3  * (x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +log(lr2 + x3 + y3) + xlogy(-1, lr1 + x3 -
                y3) + xlogy(2  * (1 - 2  * nu) * nu, lr2 + x3 + y3)))
        end

        function J1223d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * (-lr1 ^ (-1) * (x1 -
            y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) - ((-1) - 2  * nu +
            4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (
            -1) -4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) - 4  * lr2 ^ (
            -1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 2  * (x1 + (-1)
                         * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            (x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu *
            (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  *
            x3) +(x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) *
            y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 +
            (-1) * (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 + (-1)
                         * y2) * ((-3) * nu * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3)
                        ))) * G ^ (-1))
        end

        function J1223d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((-x1 + y1) * (lr1 + x3 + (-1)
                         * y3) ^ (-1) - lr1 ^ (-1) * (x1 - y1) * (x3 - y3) * (lr1 + x3 + (
            -1) * y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 - ((-1) - 2  * nu + 4  *
            nu ^ 2) * (x1 - y1) * (lr2 + x3 + y3) ^ (-1) - ((-1) - 2  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 2  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * y3 * (x3 + y3) * (2  * x3 + y3) - 2  * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (3  * nu * (
            (-3) + 2  * nu) * x3 ^ 2 + nu * ((-3) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) +2  * (2 - 9  * nu + 6  * nu ^ 2) * x3 * y3 + (3 - 9  * nu + 6  * nu ^ 2) *
            y3 ^ 2) +2  * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2)
                         * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * (
            (-2) * lr2 + 3  * x3) +(x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 + (
            -3) * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 - (lr2 - x3 - y3) * (2  * x3 + y3)))) * G ^ (
            -1))
        end

        function J2212d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x3 ^ 2  * (x3 ^ 2 + (x1 +
            (-1) * y1) ^ 2) ^ (-1) + 9  * x3 ^ 2  * (9  * x3 ^ 2 + (x1 - y1) ^ 2) ^ (-1) +
            4  * nu ^ 2  * x3 ^ 2  * (nu ^ 2  * x3 ^ 2 + (x1 - y1) ^ 2) ^ (-1) - ((-3)
            +4  * nu) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) - ((-3) + 4  * nu)
                         * lr1 ^ (-1) * (x1 - y1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) + (
            5 + 4  * nu * ((-3) + 2  * nu)) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) + (5 +
            4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 - y1) * (lr2 + x1 - y1) ^ (
            -1) * (x2 - y2) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * (
            lr1 + x2 - y2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1)
                         ^ 2  * (lr2 + x2 - y2) ^ (-1) + 4  * ((-1) + nu) * lr1 * (x2 - y2) * ((x1 +
            (-1) * y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 -
            y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 + (
            -1) * y1) ^ 2  * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (
            -1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 + (
            -4) * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         ^ (-1) * (x2 - y2) * (x3 + y3) + y3 ^ 2  * ((x1 - y1) ^ 2 + y3 ^ 2) ^ (
            -1) +9  * y3 ^ 2  * ((x1 - y1) ^ 2 + 9  * y3 ^ 2) ^ (-1) + 4  * nu ^ 2  * y3 ^ 2  * (
            (x1 - y1) ^ 2 + nu ^ 2  * y3 ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * x3 * (x2 -
            y2) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) * (x3 + y3) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) +
            nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 * (x2 - y2) * (x3 + y3)
                         ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (x2 + (
            -1) * y2) * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) ^ 2  * (x2 + (-1)
                         * y2) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 +
            (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + xlogy(4 - 4  * nu, lr1 + x2 -
                y2) + xlogy(4 - 4  * nu, lr2 + x2 - y2)))
        end

        function J2212d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (-((-3) + 4  * nu)
                         * lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) ^ 2 + (5 + 4  * nu * ((
                        - 3) + 2  * nu)) * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) ^ 2 + (
            -4) * ((-1) + nu) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1) - 4  * ((-1) +
            nu) * lr1 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr1 + x2 - y2) ^ (
            -1) -4  * ((-1) + nu) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1) - 4  * (
            (-1) + nu) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x2 - y2)
                         ^ (-1) + 4  * ((-1) + nu) * lr1 * (x1 - y1) * ((x1 - y1) ^ 2 + (x3 + (-1)
                         * y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 + (
            -1) * y3) ^ 2 - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) * (x2 -
            y2) ^ 2  * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2)
                         ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 + 4  * ((-1) + nu) * ((-1) +
            2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) *
            (x3 + y3) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) ^ 2  * y3 * (
            (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * x3 * (x1 - y1)
                         * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x3 + y3) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) +
            nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) ^ 2  * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 * (x1 - y1) * (x3 + y3)
                         ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * (x2 -
            y2) ^ 2  * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) * (x2 - y2)
                         ^ 2  * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + xlogy(3 - 4  * nu, lr1 + x1 -
                y1) + xlogy(5 + 4  * nu * ((-3) + 2  * nu), lr2 + x1 - y1)))
        end

        function J2212d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x3 * (x3 ^ 2 + (x1 + (
            -1) * y1) ^ 2) ^ (-1) * (-x1 + y1) + 9  * x3 * (9  * x3 ^ 2 + (x1 - y1)
                         ^ 2) ^ (-1) * (-x1 + y1) + 4  * nu ^ 2  * x3 * (nu ^ 2  * x3 ^ 2 + (x1 -
            y1) ^ 2) ^ (-1) * (-x1 + y1) - ((-3) + 4  * nu) * lr1 ^ (-1) * (lr1 + x1 +
            (-1) * y1) ^ (-1) * (x2 - y2) * (x3 - y3) - 4  * ((-1) + nu) *
            lr1 ^ (-1) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + (
            -4) * ((-1) + nu) * lr1 * (x1 - y1) * (x2 - y2) * ((x1 - y1)
                         ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) * (
            x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 3 + (5 + 4  * nu * ((-3)
            +2  * nu)) * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 + y3)
            -4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1)
                         * (x3 + y3) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (
            x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * (x1 + (-1)
                         * y1) * (x2 - y2) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * (
            (-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            x3 + y3) ^ 3  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) *
            lr2 * (x1 - y1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu)
                         * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (x3 + y3) ^ 3  * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)
            -2  * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1) * (
                x2 - y2) ^ (-1)) -atan(x3, x1 - y1) - 3  * atan(3  * x3,
                    x1 - y1) + 4  * nu * atan(-nu * x3, x1 - y1) + 4  * ((-1) + nu) *
            ((-1) + 2  * nu) * atan(lr2 * (-x2 + y2), (x1 - y1) * (x3 + y3)) + (4 + (
            -4) * nu) * atan(lr1 * (x3 - y3), (x1 - y1) * (x2 - y2)) + (
            4 - 4  * nu) * atan(lr2 * (x3 + y3), (x1 - y1) * (x2 - y2))))
        end

        function J2213d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x2 ^ 2  * (x2 ^ 2 + (x1 +
            (-1) * y1) ^ 2) ^ (-1) + 9  * x2 ^ 2  * (9  * x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) +
            4  * nu ^ 2  * x2 ^ 2  * (nu ^ 2  * x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) + 2  * x3 * (lr2 +
            (-1) * x1 + y1) ^ (-1) + 2  * lr2 ^ (-1) * x3 * (-x1 + y1) * (lr2 - x1 +
            y1) ^ (-1) - 2  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2 + y2 ^ 2  * ((x1 - y1) ^ 2 +
            y2 ^ 2) ^ (-1) + 9  * y2 ^ 2  * ((x1 - y1) ^ 2 + 9  * y2 ^ 2) ^ (-1) + 4  *
            nu ^ 2  * y2 ^ 2  * ((x1 - y1) ^ 2 + nu ^ 2  * y2 ^ 2) ^ (-1) - ((-3) +
            4  * nu) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) *
            lr1 ^ (-1) * (x1 - y1) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) + 2  *
            ((-1) + 2  * nu) * lr1 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 -
            y3) -2  * ((-1) + 2  * nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (
            x3 - y3) ^ 2) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) * lr1 ^ (-1)
                         * (x1 - y1) ^ 2  * (lr1 + x3 - y3) ^ (-1) - (5 + 4  * nu * ((-3) +
            2  * nu)) * (lr2 + x1 - y1) ^ (-1) * (x3 + y3) - (5 + 4  * nu * ((-3) + 2  *
            nu)) * lr2 ^ (-1) * (x1 - y1) * (lr2 + x1 - y1) ^ (-1) * (x3 + y3) + (
            -1) * (3 - 6  * nu + 4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x3 + y3)
                         ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * y3 * (2  * x3 + y3) - 2  * ((-1) + nu) * ((
                        - 1) +2  * nu) * lr2 ^ (-2) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * y3 * (2  * x3 + y3) + 2  * ((-1) + 2  * nu) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) *
            (((-1) + 2  * nu) * lr2 * (x2 - y2) ^ 2  * (x3 + y3) - ((-1) + nu) * y3 *
            (2  * x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)) -4  * lr2 ^ (-1) * (x1 + (
            -1) * y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (x2 ^ 4  * y3 + 4  * x2 ^ 2  * x3 ^ 2  * y3 + 2  *
            x3 ^ 4  * y3 + x2 ^ 2  * y1 ^ 2  * y3 + 2  * x3 ^ 2  * y1 ^ 2  * y3 - 4  * x2 ^ 3  * y2 *
            y3 - 8  * x2 * x3 ^ 2  * y2 * y3 - 2  * x2 * y1 ^ 2  * y2 * y3 + 6  * x2 ^ 2  *
            y2 ^ 2  * y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 + y1 ^ 2  * y2 ^ 2  * y3 - 4  * x2 * y2 ^ 3  *
            y3 + y2 ^ 4  * y3 + 6  * x2 ^ 2  * x3 * y3 ^ 2 + 7  * x3 ^ 3  * y3 ^ 2 + 3  * x3 * y1 ^ 2  *
            y3 ^ 2 - 12  * x2 * x3 * y2 * y3 ^ 2 + 6  * x3 * y2 ^ 2  * y3 ^ 2 + 2  * x2 ^ 2  *
            y3 ^ 3 + 9  * x3 ^ 2  * y3 ^ 3 + y1 ^ 2  * y3 ^ 3 - 4  * x2 * y2 * y3 ^ 3 + 2  * y2 ^ 2  *
            y3 ^ 3 + 5  * x3 * y3 ^ 4 + y3 ^ 5 - 3  * nu * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) + 2  * nu ^ 2  *
            (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) +x1 ^ 2  * y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 +
            y3)) -2  * x1 * y1 * y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +(
            -2) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) *
            ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (x2 ^ 4  * y3 - 2  * lr2 * x2 ^ 2  * x3 * y3 + 4  *
            x2 ^ 2  * x3 ^ 2  * y3 - 2  * lr2 * x3 ^ 3  * y3 + 2  * x3 ^ 4  * y3 + x2 ^ 2  * y1 ^ 2  *
            y3 + 2  * x3 ^ 2  * y1 ^ 2  * y3 - 4  * x2 ^ 3  * y2 * y3 + 4  * lr2 * x2 * x3 * y2 * y3 + (
            -8) * x2 * x3 ^ 2  * y2 * y3 - 2  * x2 * y1 ^ 2  * y2 * y3 + 6  * x2 ^ 2  * y2 ^ 2  *
            y3 - 2  * lr2 * x3 * y2 ^ 2  * y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 + y1 ^ 2  * y2 ^ 2  * y3 + (
            -4) * x2 * y2 ^ 3  * y3 + y2 ^ 4  * y3 - lr2 * x2 ^ 2  * y3 ^ 2 + 6  * x2 ^ 2  * x3 *
            y3 ^ 2 - 5  * lr2 * x3 ^ 2  * y3 ^ 2 + 7  * x3 ^ 3  * y3 ^ 2 + 3  * x3 * y1 ^ 2  * y3 ^ 2 +
            2  * lr2 * x2 * y2 * y3 ^ 2 - 12  * x2 * x3 * y2 * y3 ^ 2 - lr2 * y2 ^ 2  *
            y3 ^ 2 + 6  * x3 * y2 ^ 2  * y3 ^ 2 + 2  * x2 ^ 2  * y3 ^ 3 - 4  * lr2 * x3 * y3 ^ 3 + 9  *
            x3 ^ 2  * y3 ^ 3 + y1 ^ 2  * y3 ^ 3 - 4  * x2 * y2 * y3 ^ 3 + 2  * y2 ^ 2  * y3 ^ 3 + (-1)
                         * lr2 * y3 ^ 4 + 5  * x3 * y3 ^ 4 + y3 ^ 5 - 3  * nu * (x3 * (x3 ^ 2 + y1 ^ 2 + (x2 + (
            -1) * y2) ^ 2) +((-2) * lr2 * x3 + 3  * x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) * y3 + (
            -1) * (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) +
            2  * nu ^ 2  * (x3 * (x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) + ((-2) * lr2 * x3 + 3  *
            x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + 2  * x1 * y1 * (3  * nu * (x3 + y3) *
            ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) - 2  * nu ^ 2  * (x3 + y3) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) -y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +
            x1 ^ 2  * ((-3) * nu * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + 2  *
            nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + y3 * ((x2 - y2)
                         ^ 2 + (x3 + y3) * (2  * x3 + y3)))) +2  * lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (-x3 *
            y1 ^ 2  * y2 ^ 2 + 2  * x3 ^ 4  * y3 + 6  * x3 ^ 2  * y1 ^ 2  * y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 +
            2  * y1 ^ 2  * y2 ^ 2  * y3 + y2 ^ 4  * y3 + 7  * x3 ^ 3  * y3 ^ 2 + 9  * x3 * y1 ^ 2  * y3 ^ 2 +
            6  * x3 * y2 ^ 2  * y3 ^ 2 + 9  * x3 ^ 2  * y3 ^ 3 + 3  * y1 ^ 2  * y3 ^ 3 + 2  * y2 ^ 2  *
            y3 ^ 3 + 5  * x3 * y3 ^ 4 + y3 ^ 5 + x2 ^ 4  * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (
            x3 + y3)) -4  * x2 ^ 3  * y2 * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) +
            2  * x2 * y2 * (x3 * y1 ^ 2 - 4  * x3 ^ 2  * y3 - 2  * y1 ^ 2  * y3 - 2  *
            y2 ^ 2  * y3 - 6  * x3 * y3 ^ 2 - 2  * y3 ^ 3 - 2  * nu ^ 2  * (x3 + y3) * (
            y1 ^ 2 + 2  * y2 ^ 2 + 2  * (x3 + y3) ^ 2) +nu * (x3 + y3) * (5  * y1 ^ 2 + 6  * y2 ^ 2 + 6  * (
            x3 + y3) ^ 2)) +x2 ^ 2  * (-x3 * y1 ^ 2 + 4  * x3 ^ 2  * y3 + 2  * y1 ^ 2  * y3 + 6  *
            y2 ^ 2  * y3 + 6  * x3 * y3 ^ 2 + 2  * y3 ^ 3 + 2  * nu ^ 2  * (x3 + y3) * (y1 ^ 2 + 6  *
            y2 ^ 2 + 2  * (x3 + y3) ^ 2) -nu * (x3 + y3) * (5  * y1 ^ 2 + 18  * y2 ^ 2 + 6  * (
            x3 + y3) ^ 2)) +x1 ^ 2  * (-x2 ^ 2  * x3 + 2  * x2 * x3 * y2 - x3 * y2 ^ 2 +
            2  * x2 ^ 2  * y3 + 6  * x3 ^ 2  * y3 - 4  * x2 * y2 * y3 + 2  * y2 ^ 2  * y3 + 9  * x3 *
            y3 ^ 2 + 3  * y3 ^ 3 + 2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + 3  * (x3 + y3) ^ 2)
            -nu * (x3 + y3) * (5  * (x2 - y2) ^ 2 + 9  * (x3 + y3) ^ 2)) +2  * x1 *
            y1 * (x2 ^ 2  * x3 - 2  * x2 * x3 * y2 + x3 * y2 ^ 2 - 2  * x2 ^ 2  * y3 - 6  *
            x3 ^ 2  * y3 + 4  * x2 * y2 * y3 - 2  * y2 ^ 2  * y3 - 9  * x3 * y3 ^ 2 - 3  *
            y3 ^ 3 - 2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + 3  * (x3 + y3) ^ 2) + nu * (
            x3 + y3) * (5  * (x2 - y2) ^ 2 + 9  * (x3 + y3) ^ 2)) -nu * (x3 + y3) * (
            3  * x3 ^ 4 + 12  * x3 ^ 3  * y3 + 3  * (y2 ^ 2 + y3 ^ 2) ^ 2 + 3  * x3 ^ 2  * (3  * y1 ^ 2 + 2  *
            y2 ^ 2 + 6  * y3 ^ 2) +y1 ^ 2  * (5  * y2 ^ 2 + 9  * y3 ^ 2) + 6  * x3 * y3 * (3  * y1 ^ 2 +
            2  * (y2 ^ 2 + y3 ^ 2))) +2  * nu ^ 2  * (x3 + y3) * (x3 ^ 4 + 4  * x3 ^ 3  * y3 + (y2 ^ 2 +
            y3 ^ 2) ^ 2 + y1 ^ 2  * (y2 ^ 2 + 3  * y3 ^ 2) + x3 ^ 2  * (3  * y1 ^ 2 + 2  * y2 ^ 2 + 6  *
            y3 ^ 2) +x3 * (6  * y1 ^ 2  * y3 + 4  * y3 * (y2 ^ 2 + y3 ^ 2)))) +xlogy(3 - 4  * nu,
                lr1 + x3 - y3) + xlogy((-3) + 6  * nu - 4  * nu ^ 2, lr2 + x3 + y3)))
        end

        function J2213d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x2 * (x2 ^ 2 + (x1 + (
            -1) * y1) ^ 2) ^ (-1) * (-x1 + y1) + 9  * x2 * (9  * x2 ^ 2 + (x1 - y1)
                         ^ 2) ^ (-1) * (-x1 + y1) + 4  * nu ^ 2  * x2 * (nu ^ 2  * x2 ^ 2 + (x1 -
            y1) ^ 2) ^ (-1) * (-x1 + y1) + 2  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 + (-1)
                         * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) +
            2  * lr2 ^ (-1) * x3 * (lr2 - x1 + y1) ^ (-1) * (-x2 + y2) - ((-3)
            +4  * nu) * lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 + (-1)
                         * y3) -2  * ((-1) + 2  * nu) * lr1 * (x1 - y1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * ((x2 - y2) ^ 2 + (x3 + (-1)
                         * y3) ^ 2) ^ (-1) * (x3 - y3) - 2  * ((-1) + 2  * nu) * lr1 ^ (-1) * (x1 + (
            -1) * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2)
                         ^ 3  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) + (-1)
                         * ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr1 + x3 + (
            -1) * y3) ^ (-1) - (5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (lr2 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) * (x3 + y3) - (3 - 6  * nu + 4  * nu ^ 2)
                         * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + 4  * ((
                        - 1) +nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3) - 2  * ((-1) + nu) * ((
                        - 1) +2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) + 4  * ((-1) + nu) * ((-1) +
            2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) *
            (x2 - y2) * y3 * (2  * x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)
            -2  * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)
                         * (((-1) + 2  * nu) * lr2 * (x3 + y3) + 2  * ((-1) + nu) * y3 * (2  * x3 + y3)) + 2  *
            lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x2 - y2) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (-
            x2 ^ 2  * x3 + 2  * x2 * x3 * y2 - x3 * y2 ^ 2 + 2  * x1 ^ 2  * y3 + 3  * x2 ^ 2  * y3 +
            8  * x3 ^ 2  * y3 - 4  * x1 * y1 * y3 + 2  * y1 ^ 2  * y3 - 6  * x2 * y2 * y3 + 3  *
            y2 ^ 2  * y3 + 12  * x3 * y3 ^ 2 + 4  * y3 ^ 3 + 4  * nu ^ 2  * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3) ^ 2) -2  * nu * (x3 + y3) * (3  * (x1 + (
            -1) * y1) ^ 2 + 4  * (x2 - y2) ^ 2 + 6  * (x3 + y3) ^ 2)) -4  * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) * (x2 ^ 4  * y3 + 4  * x2 ^ 2  *
            x3 ^ 2  * y3 + 2  * x3 ^ 4  * y3 + x2 ^ 2  * y1 ^ 2  * y3 + 2  * x3 ^ 2  * y1 ^ 2  * y3 - 4  *
            x2 ^ 3  * y2 * y3 - 8  * x2 * x3 ^ 2  * y2 * y3 - 2  * x2 * y1 ^ 2  * y2 * y3 + 6  *
            x2 ^ 2  * y2 ^ 2  * y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 + y1 ^ 2  * y2 ^ 2  * y3 - 4  * x2 *
            y2 ^ 3  * y3 + y2 ^ 4  * y3 + 6  * x2 ^ 2  * x3 * y3 ^ 2 + 7  * x3 ^ 3  * y3 ^ 2 + 3  * x3 *
            y1 ^ 2  * y3 ^ 2 - 12  * x2 * x3 * y2 * y3 ^ 2 + 6  * x3 * y2 ^ 2  * y3 ^ 2 + 2  *
            x2 ^ 2  * y3 ^ 3 + 9  * x3 ^ 2  * y3 ^ 3 + y1 ^ 2  * y3 ^ 3 - 4  * x2 * y2 * y3 ^ 3 + 2  *
            y2 ^ 2  * y3 ^ 3 + 5  * x3 * y3 ^ 4 + y3 ^ 5 - 3  * nu * (x3 + y3) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +
            2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +x1 ^ 2  * y3 * ((x2 - y2) ^ 2 + (x3 +
            y3) * (2  * x3 + y3)) -2  * x1 * y1 * y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  *
            x3 + y3))) -4  * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-2) * (x2 - y2) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) * (x2 ^ 4  * y3 + 4  * x2 ^ 2  * x3 ^ 2  * y3 + 2  * x3 ^ 4  * y3 + x2 ^ 2  * y1 ^ 2  * y3 +
            2  * x3 ^ 2  * y1 ^ 2  * y3 - 4  * x2 ^ 3  * y2 * y3 - 8  * x2 * x3 ^ 2  * y2 * y3 + (
            -2) * x2 * y1 ^ 2  * y2 * y3 + 6  * x2 ^ 2  * y2 ^ 2  * y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 +
            y1 ^ 2  * y2 ^ 2  * y3 - 4  * x2 * y2 ^ 3  * y3 + y2 ^ 4  * y3 + 6  * x2 ^ 2  * x3 *
            y3 ^ 2 + 7  * x3 ^ 3  * y3 ^ 2 + 3  * x3 * y1 ^ 2  * y3 ^ 2 - 12  * x2 * x3 * y2 *
            y3 ^ 2 + 6  * x3 * y2 ^ 2  * y3 ^ 2 + 2  * x2 ^ 2  * y3 ^ 3 + 9  * x3 ^ 2  * y3 ^ 3 + y1 ^ 2  *
            y3 ^ 3 - 4  * x2 * y2 * y3 ^ 3 + 2  * y2 ^ 2  * y3 ^ 3 + 5  * x3 * y3 ^ 4 + y3 ^ 5
            -3  * nu * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) + x1 ^ 2  *
            y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3)) - 2  * x1 * y1 * y3 * ((x2 +
            (-1) * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) -2  * (x1 - y1) * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) * (x2 ^ 4  * y3 - 2  * lr2 * x2 ^ 2  * x3 * y3 + 4  * x2 ^ 2  * x3 ^ 2  *
            y3 - 2  * lr2 * x3 ^ 3  * y3 + 2  * x3 ^ 4  * y3 + x2 ^ 2  * y1 ^ 2  * y3 + 2  * x3 ^ 2  *
            y1 ^ 2  * y3 - 4  * x2 ^ 3  * y2 * y3 + 4  * lr2 * x2 * x3 * y2 * y3 - 8  * x2 *
            x3 ^ 2  * y2 * y3 - 2  * x2 * y1 ^ 2  * y2 * y3 + 6  * x2 ^ 2  * y2 ^ 2  * y3 - 2  *
            lr2 * x3 * y2 ^ 2  * y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 + y1 ^ 2  * y2 ^ 2  * y3 - 4  * x2 *
            y2 ^ 3  * y3 + y2 ^ 4  * y3 - lr2 * x2 ^ 2  * y3 ^ 2 + 6  * x2 ^ 2  * x3 * y3 ^ 2 + (
            -5) * lr2 * x3 ^ 2  * y3 ^ 2 + 7  * x3 ^ 3  * y3 ^ 2 + 3  * x3 * y1 ^ 2  * y3 ^ 2 + 2  * lr2 *
            x2 * y2 * y3 ^ 2 - 12  * x2 * x3 * y2 * y3 ^ 2 - lr2 * y2 ^ 2  * y3 ^ 2 + 6  *
            x3 * y2 ^ 2  * y3 ^ 2 + 2  * x2 ^ 2  * y3 ^ 3 - 4  * lr2 * x3 * y3 ^ 3 + 9  * x3 ^ 2  *
            y3 ^ 3 + y1 ^ 2  * y3 ^ 3 - 4  * x2 * y2 * y3 ^ 3 + 2  * y2 ^ 2  * y3 ^ 3 - lr2 *
            y3 ^ 4 + 5  * x3 * y3 ^ 4 + y3 ^ 5 - 3  * nu * (x3 * (x3 ^ 2 + y1 ^ 2 + (x2 - y2)
                         ^ 2) +((-2) * lr2 * x3 + 3  * x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 +
            (-3) * x3) * y3 ^ 2 + y3 ^ 3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + 2  * nu ^ 2  * (
            x3 * (x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) + ((-2) * lr2 * x3 + 3  * x3 ^ 2 + y1 ^ 2 + (
            x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) +2  * x1 * y1 * (3  * nu * (x3 + y3) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) -2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) -y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +x1 ^ 2  * ((
                        - 3) * nu * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + 2  * nu ^ 2  * (x3 + y3)
                         * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (
            2  * x3 + y3)))) +2  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1) ^ (-1)
                         * (x2 - y2)) + atan(x1 - y1, (-1) * x2) - 3  * atan(3  * x2, x1 + (
                -1) * y1) +4  * nu * atan(-nu * x2, x1 - y1) + ((-2) + 4  * nu) *
            atan(lr1 * (-x2 + y2), (x1 - y1) * (x3 - y3)) + 2  * (1 - 2  *
            nu) ^ 2  * atan(lr2 * (-x2 + y2), (x1 - y1) * (x3 + y3))))
        end

        function J2213d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (2  * ((-1) + 2  * nu) *
            lr1 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 +
            (-1) * y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) - ((
                        - 3) +4  * nu) * lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) ^ 2 + (
            -2) * ((-1) + 2  * nu) * lr1 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 + (-1)
                         * y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 - ((-3) + 4  * nu) * (x1 -
            y1) * (lr1 + x3 - y3) ^ (-1) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 + (-1)
                         * y1) * (x3 - y3) * (lr1 + x3 - y3) ^ (-1) - 2  * lr2 ^ (-1) * x3 * (
            lr2 - x1 + y1) ^ (-1) * (x3 + y3) - (5 + 4  * nu * ((-3) + 2  * nu)) *
            lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x3 + y3) ^ 2 - (3 - 6  * nu + 4  *
            nu ^ 2) * (x1 - y1) * (lr2 + x3 + y3) ^ (-1) - (3 - 6  * nu + 4  *
            nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) - 2  * (
            (-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 * (x3 + y3) * (2  * x3 + y3) + 4  * ((-1) + nu)
                         * ((-1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         ^ (-1) * y3 * (x3 + y3) * (2  * x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) -4  * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) * (x2 ^ 4  *
            y3 + 4  * x2 ^ 2  * x3 ^ 2  * y3 + 2  * x3 ^ 4  * y3 + x2 ^ 2  * y1 ^ 2  * y3 + 2  * x3 ^ 2  *
            y1 ^ 2  * y3 - 4  * x2 ^ 3  * y2 * y3 - 8  * x2 * x3 ^ 2  * y2 * y3 - 2  * x2 *
            y1 ^ 2  * y2 * y3 + 6  * x2 ^ 2  * y2 ^ 2  * y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 + y1 ^ 2  *
            y2 ^ 2  * y3 - 4  * x2 * y2 ^ 3  * y3 + y2 ^ 4  * y3 + 6  * x2 ^ 2  * x3 * y3 ^ 2 + 7  *
            x3 ^ 3  * y3 ^ 2 + 3  * x3 * y1 ^ 2  * y3 ^ 2 - 12  * x2 * x3 * y2 * y3 ^ 2 + 6  * x3 *
            y2 ^ 2  * y3 ^ 2 + 2  * x2 ^ 2  * y3 ^ 3 + 9  * x3 ^ 2  * y3 ^ 3 + y1 ^ 2  * y3 ^ 3 - 4  *
            x2 * y2 * y3 ^ 3 + 2  * y2 ^ 2  * y3 ^ 3 + 5  * x3 * y3 ^ 4 + y3 ^ 5 - 3  * nu * (x3 + y3)
                         * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) *
            ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) + x1 ^ 2  * y3 * ((x2 + (
            -1) * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3)) -2  * x1 * y1 * y3 * ((x2 - y2)
                         ^ 2 + (x3 + y3) * (2  * x3 + y3))) +2  * ((-1) + 2  * nu) * (x1 - y1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-1) * (((-1) + 2  * nu) * lr2 * (x2 - y2) ^ 2 - 2  * ((-1) + nu) * y3 * (
            (x2 - y2) ^ 2 + (x3 + y3) * (3  * x3 + 2  * y3))) +2  * lr2 ^ (-1) * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 +
            (x3 + y3) ^ 2) ^ (-1) * (nu * ((-3) + 2  * nu) * x2 ^ 4 + 4  * (3 - 2  * nu) * nu *
            x2 ^ 3  * y2 - x3 ^ 2  * y2 ^ 2 + 4  * x1 ^ 2  * x3 * y3 + 8  * x3 ^ 3  * y3 - 8  *
            x1 * x3 * y1 * y3 + 4  * x3 * y1 ^ 2  * y3 + 6  * x3 * y2 ^ 2  * y3 + 3  * x1 ^ 2  * y3 ^ 2 +
            21  * x3 ^ 2  * y3 ^ 2 - 6  * x1 * y1 * y3 ^ 2 + 3  * y1 ^ 2  * y3 ^ 2 + 5  * y2 ^ 2  *
            y3 ^ 2 + 18  * x3 * y3 ^ 3 + 5  * y3 ^ 4 + 2  * x2 * y2 * ((1 + 2  * (7 - 4  * nu) * nu) *
            x3 ^ 2 + (3 - 2  * nu) * nu * ((x1 - y1) ^ 2 + 2  * y2 ^ 2) - 2  * (3 + 2  *
            nu * ((-7) + 4  * nu)) * x3 * y3 + ((-5) + 2  * (7 - 4  * nu) * nu) * y3 ^ 2) +
            x2 ^ 2  * (((-1) + 2  * nu * ((-7) + 4  * nu)) * x3 ^ 2 + nu * ((-3) + 2  * nu) * ((x1 + (
            -1) * y1) ^ 2 + 6  * y2 ^ 2) +2  * (3 + 2  * nu * ((-7) + 4  * nu)) * x3 * y3 + (5 + 2  *
            nu * ((-7) + 4  * nu)) * y3 ^ 2) -nu * (y2 ^ 2 + 3  * (x3 + y3) ^ 2) * (3  * (
            x1 - y1) ^ 2 + 3  * y2 ^ 2 + 5  * (x3 + y3) ^ 2) +2  * nu ^ 2  * (5  * x3 ^ 4 +
            y1 ^ 2  * y2 ^ 2 + y2 ^ 4 + 20  * x3 ^ 3  * y3 + 3  * y1 ^ 2  * y3 ^ 2 + 4  * y2 ^ 2  * y3 ^ 2 +
            5  * y3 ^ 4 + x3 ^ 2  * (3  * y1 ^ 2 + 4  * y2 ^ 2 + 30  * y3 ^ 2) + x3 * (6  * y1 ^ 2  * y3 +
            8  * y2 ^ 2  * y3 + 20  * y3 ^ 3) +x1 ^ 2  * (y2 ^ 2 + 3  * (x3 + y3) ^ 2) - 2  * x1 *
            y1 * (y2 ^ 2 + 3  * (x3 + y3) ^ 2))) -2  * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2)
                         * (x2 ^ 4  * y3 - 2  * lr2 * x2 ^ 2  * x3 * y3 + 4  * x2 ^ 2  * x3 ^ 2  * y3 - 2  *
            lr2 * x3 ^ 3  * y3 + 2  * x3 ^ 4  * y3 + x2 ^ 2  * y1 ^ 2  * y3 + 2  * x3 ^ 2  * y1 ^ 2  * y3 + (
            -4) * x2 ^ 3  * y2 * y3 + 4  * lr2 * x2 * x3 * y2 * y3 - 8  * x2 * x3 ^ 2  * y2 * y3 + (
            -2) * x2 * y1 ^ 2  * y2 * y3 + 6  * x2 ^ 2  * y2 ^ 2  * y3 - 2  * lr2 * x3 * y2 ^ 2  *
            y3 + 4  * x3 ^ 2  * y2 ^ 2  * y3 + y1 ^ 2  * y2 ^ 2  * y3 - 4  * x2 * y2 ^ 3  * y3 +
            y2 ^ 4  * y3 - lr2 * x2 ^ 2  * y3 ^ 2 + 6  * x2 ^ 2  * x3 * y3 ^ 2 - 5  * lr2 *
            x3 ^ 2  * y3 ^ 2 + 7  * x3 ^ 3  * y3 ^ 2 + 3  * x3 * y1 ^ 2  * y3 ^ 2 + 2  * lr2 * x2 * y2 *
            y3 ^ 2 - 12  * x2 * x3 * y2 * y3 ^ 2 - lr2 * y2 ^ 2  * y3 ^ 2 + 6  * x3 *
            y2 ^ 2  * y3 ^ 2 + 2  * x2 ^ 2  * y3 ^ 3 - 4  * lr2 * x3 * y3 ^ 3 + 9  * x3 ^ 2  * y3 ^ 3 +
            y1 ^ 2  * y3 ^ 3 - 4  * x2 * y2 * y3 ^ 3 + 2  * y2 ^ 2  * y3 ^ 3 - lr2 * y3 ^ 4 +
            5  * x3 * y3 ^ 4 + y3 ^ 5 - 3  * nu * (x3 * (x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) + (
            (-2) * lr2 * x3 + 3  * x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  *
            x3) * y3 ^ 2 + y3 ^ 3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + 2  * nu ^ 2  * (x3 * (
            x3 ^ 2 + y1 ^ 2 + (x2 - y2) ^ 2) +((-2) * lr2 * x3 + 3  * x3 ^ 2 + y1 ^ 2 + (x2 + (
            -1) * y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) +2  * x1 * y1 * (3  * nu * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) -2  * nu ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) -
            y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +x1 ^ 2  * ((-3) * nu * (x3 +
            y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) + 2  * nu ^ 2  * (x3 + y3) * ((x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x2 - y2) ^ 2 + (x3 + y3) * (2  * x3 + y3)))) +
            xlogy((-2), lr2 - x1 + y1) + xlogy(3 - 4  * nu, lr1 + x1 - y1) + xlogy(
                (-5) + 4  * (3 - 2  * nu) * nu, lr2 + x1 - y1)))
        end

        function J2223d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x1 * (x1 ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (-x2 + y2) + 9  * x1 * (9  * x1 ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (-x2 + y2) + 4  * nu ^ 2  * x1 * (nu ^ 2  * x1 ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (-x2 + y2) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 -
            y1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * lr1 * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) + (
            -4) * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 3  * ((x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 -
            y3) ^ 2) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 + (
            -1) * y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) + 4  * ((-1) + nu) *
            lr2 ^ (-1) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) - (3 +
            (-6) * nu + 4  * nu ^ 2) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 +
            x3 + y3) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 - y1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-2) * (x2 - y2) * y3 * (2  * x3 + y3)
            +2  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3)
            -4  * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) * (y3 - 3  * nu * (x3 + y3) + 2  * nu ^ 2  * (x3 + y3)) + 4  * ((
                        - 1) +nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) ^ 3  * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) + 2  * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 * (x3 ^ 2 + (x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +y3 * ((
                x1- y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 - y3) * (
            2  * x3 + y3))) +4  * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-2) * (x2 - y2) * ((-3) * nu * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +atan(x2 - y2, (-1) * x1) - 3  *
            atan(3  * x1, x2 - y2) + 4  * nu * atan(-nu * x1, x2 - y2) + (
            4 - 4  * nu) * atan(lr1 * (x1 - y1), (x2 - y2) * (x3 - y3))
            +4  * ((-1) + nu) * atan(lr2 * (x1 - y1), (x2 - y2) * (x3 + y3))))
        end

        function J2223d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (x1 ^ 2  * (x1 ^ 2 + (x2 +
            (-1) * y2) ^ 2) ^ (-1) + 9  * x1 ^ 2  * (9  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) +
            4  * nu ^ 2  * x1 ^ 2  * (nu ^ 2  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + y1 ^ 2  * (
            y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 9  * y1 ^ 2  * (9  * y1 ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) + 4  * nu ^ 2  * y1 ^ 2  * (nu ^ 2  * y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) +
            (-4) * ((-1) + nu) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) - 4  * ((-1)
            +nu) * lr1 ^ (-1) * (x2 - y2) * (lr1 + x2 - y2) ^ (-1) * (x3 -
            y3) +4  * ((-1) + nu) * lr1 * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 + (
            -1) * y3) -4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x1 - y1)
                         ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) - ((-3) + 4  * nu) *
            lr1 ^ (-1) * (x2 - y2) ^ 2  * (lr1 + x3 - y3) ^ (-1) + 4  * ((-1) + nu) * (
            lr2 + x2 - y2) ^ (-1) * (x3 + y3) + 4  * ((-1) + nu) * lr2 ^ (-1) * (x2 -
            y2) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) - (3 - 6  * nu + 4  * nu ^ 2) *
            lr2 ^ (-1) * (x2 - y2) ^ 2  * (lr2 + x3 + y3) ^ (-1) + 2  * ((-1) + nu) * ((-1) +
            2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * y3 * (2  * x3 + y3) +
            (-4) * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         ^ (-2) * (x2 - y2) ^ 2  * y3 * (2  * x3 + y3) + 2  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) ^ 2  * y3 * (2  * x3 + y3) - 4  * ((-1) + nu) * lr2 * (x1 - y1) ^ 2  * ((
                x1- y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((
                    x1- y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3)
                         * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 *
            (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +4  * lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-2) * (x2 - y2) ^ 2  * ((-3) * nu * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +2  * nu ^ 2  * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) +y3 * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +lr2 ^ (-1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (6  * nu * (x3 + y3) * ((x1 - y1) ^ 2 + 3  * (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) -4  * nu ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + 3  * (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) -2  * y3 * ((x1 - y1) ^ 2 + 3  * (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) * (2  * x3 + y3))) +xlogy(3 - 4  * nu, lr1 + x3 - y3)
            +xlogy((-3) + 6  * nu - 4  * nu ^ 2, lr2 + x3 + y3)))
        end

        function J2223d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * lr1 *
            (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) - 4  * ((-1) +
            nu) * lr1 ^ (-1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) ^ 2 - 4  * ((
                        - 1) +nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) ^ 2 - ((-3) + 4  * nu) * (x2 - y2) * (lr1 +
            x3 - y3) ^ (-1) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x2 - y2) * (
            x3 - y3) * (lr1 + x3 - y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (
            (x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * y3 + 4  * ((
                        - 1) +nu) * lr2 ^ (-1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) ^ 2 - (3 + (
            -6) * nu + 4  * nu ^ 2) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) - (3 - 6  *
            nu + 4  * nu ^ 2) * lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) +
            2  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-2) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) * y3 * (x3 + y3) * (2  * x3 + y3) - 2  * lr2 ^ (
            -1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            3  * nu * ((-3) + 2  * nu) * x3 ^ 2 + nu * ((-3) + 2  * nu) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2) +2  * (2 - 9  * nu + 6  * nu ^ 2) * x3 * y3 + (3 - 9  * nu + 6  *
            nu ^ 2) * y3 ^ 2) -4  * ((-1) + nu) * lr2 * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * ((
                x1- y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((-3) * nu * (x3 * (x3 ^ 2 + (x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) +(x3 * ((-2) * lr2 + 3  * x3) + (x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 + y3 ^ 3) +2  * nu ^ 2  * (x3 *
            (x3 ^ 2 + (x1 - y1) ^ 2 + (x2 - y2) ^ 2) + (x3 * ((-2) * lr2 + 3  * x3) + (
            x1 - y1) ^ 2 + (x2 - y2) ^ 2) * y3 - (lr2 - 3  * x3) * y3 ^ 2 +
            y3 ^ 3) +y3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 - (lr2 - x3 + (
            -1) * y3) * (2  * x3 + y3))) +xlogy(4 - 4  * nu, lr1 + x2 - y2) + xlogy(4  * (
                (-1) + nu), lr2 + x2 - y2)))
        end

        function J3212d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * ((-4) * ((-1) + nu) * (
            (-1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) ^ 2 + (lr1 + x1 - y1) ^ (-1) * (x3 - y3) + lr1 ^ (-1) * (x1 -
            y1) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x3 + y3) ^ (-1) + 4  * ((-1) + nu)
                         * ((-1) + 2  * nu) * lr2 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (
            x2 - y2) ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)
            -4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * lr2 ^ (-1) * x3 * y3 * (x3 + y3) * ((
                x2- y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * x3 * (x1 - y1) ^ 2  * y3 * (
            x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - lr2 ^ (-1) * (x1 - y1)
                         * (lr2 + x1 - y1) ^ (-1) * (x3 + 7  * y3 + 8  * nu ^ 2  * (x3 + y3) - 8  * nu * (
            x3 + 2  * y3)) +(lr2 + x1 - y1) ^ (-1) * (-x3 - 7  * y3 - 8  *
            nu ^ 2  * (x3 + y3) + 8  * nu * (x3 + 2  * y3)) +xlogy((-4) * ((-1) + nu) * ((-1) + 2  *
                nu), lr2 + x3 + y3)))
        end

        function J3212d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * ((
                        - 1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) + lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2)
                         * (x3 - y3) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 + (-1)
                         * y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + 4  * lr2 ^ (-1) * x3 * (x1 -
            y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -2) -4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1)
                         * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) ^ 3  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * x3 * (
            x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -3 / 2) -lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 +
            7  * y3 + 8  * nu ^ 2  * (x3 + y3) - 8  * nu * (x3 + 2  * y3)) +4  * ((-1) + nu) * ((-1) +
            2  * nu) * atan((x1 - y1) ^ (-1) * (x2 - y2)) + 4  * ((-1) + nu) * ((
                        - 1) +2  * nu) * atan(lr2 * (-x2 + y2), (x1 - y1) * (x3 + y3))))
        end

        function J3212d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x3 - y3) ^ 2 - 4  * ((-1) + nu) * ((-1) + 2  * nu) * (
            x1 - y1) * (lr2 + x3 + y3) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) *
            lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + 4  * lr2 ^ (-1) *
            x3 * (x1 - y1) * y3 * (x3 + y3) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-2) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) ^ 2  * (x3 + y3) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  *
            lr2 ^ (-1) * (x1 - y1) * y3 * (2  * x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 2  * x3 * (x1 - y1) * y3 * (x3 + y3) ^ 2  * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) - lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x3 + y3) * (x3 +
            7  * y3 + 8  * nu ^ 2  * (x3 + y3) - 8  * nu * (x3 + 2  * y3)) +log(lr1 + x1 - y1) +
            xlogy((-1) - 8  * ((-1) + nu) * nu, lr2 + x1 - y1)))
        end

        function J3213d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((lr1 + x1 - y1) ^ (-1) * (
            x2 - y2) +lr1 ^ (-1) * (x1 - y1) * (lr1 + x1 - y1) ^ (-1) * (x2 +
            (-1) * y2) -(1 + 8  * ((-1) + nu) * nu) * (lr2 + x1 - y1) ^ (-1) * (x2 +
            (-1) * y2) -(1 + 8  * ((-1) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) * (
            lr2 + x1 - y1) ^ (-1) * (x2 - y2) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (
            (x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) +
            2  * lr2 ^ (-1) * x3 * (x2 - y2) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-1) - 4  * ((-1) + 2  * nu) * lr2 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         ^ (-1) * (x2 - y2) * (x3 + y3) * (nu * x3 + ((-1) + nu) * y3) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 -
            y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) *
            (x3 + y3) * (nu * x3 + ((-1) + nu) * y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) -2  * x3 * (x1 - y1) ^ 2  * (x2 - y2) * y3 * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2)) * G ^ (-1))
        end

        function J3213d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) ^ 2 - (1 + 8  * ((-1) + nu) * nu) * lr2 ^ (
            -1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) ^ 2 - 4  * ((-1) + nu) * ((
                        - 1) +2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x3 + y3) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) ^ 2  *
            y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * x3 * (x1 + (-1)
                         * y1) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + 2  * nu) *
            lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 +
            y3) * (nu * x3 + ((-1) + nu) * y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            4  * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * (nu * x3 + ((-1) + nu) *
            y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) *
            (x2 - y2) ^ 2  * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + log(lr1 + x1 -
                y1) + xlogy((-1) - 8  * ((-1) + nu) * nu, lr2 + x1 - y1)))
        end

        function J3213d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) * (x3 - y3) - (1 + 8  * ((-1) + nu)
                         * nu) * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 + y3) + (
            -4) * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((
                x2- y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * (x1 - y1) * (x2 +
            (-1) * y2) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + 2  *
            nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x2 - y2) * (nu * x3 + ((-1) + nu) * y3) * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 4  * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * (nu *
            x3 + ((-1) + nu) * y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (
            x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -3 / 2) +4  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1) * (x2 -
                y2) ^ (-1)) +4  * nu * ((-1) + 2  * nu) * atan(lr2 * (x2 - y2), (x1 -
                    y1) * (x3 + y3))))
        end

        function J3223d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x1 - y1) + (
            -1) * (1 + 8  * ((-1) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) + 2  * lr2 ^ (-1) * (
            x1 - y1) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) +2  * x3 * (x1 - y1) * y3 * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 2  * ((-3) + 4  * nu) * x3 * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2)) * G ^ (-1))
        end

        function J3223d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x2 - y2) + (
            -1) * (1 + 8  * ((-1) + nu) * nu) * lr2 ^ (-1) * (x2 - y2) + 2  * lr2 ^ (-1) * (
            x2 - y2) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) +2  * x3 * (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 2  * ((-3) + 4  * nu) * x3 * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2)) * G ^ (-1))
        end

        function J3223d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x3 + (-1)
                         * y3) +lr2 ^ (-1) * (-x3 - 3  * y3 + 8  * nu * (x3 + y3) - 8  * nu ^ 2  * (
            x3 + y3)) +2  * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) +2  * lr2 ^ (-1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  *
            y3 - 6  * nu * (x3 + y3) + 4  * nu ^ 2  * (x3 + y3)) +2  * x3 * y3 * (x3 + y3) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 2  * ((-3) + 4  *
            nu) * x3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) ^ 2  * ((
                x1- y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + 2  * ((-3) + 4  *
            nu) * lr2 ^ (-1) * x3 * (1 - (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)) ^ (-1) + 2  * ((-3) + 4  * nu) * real(acoth(complex(lr2 ^ (-1)
                         * (x3 + y3)))) + xlogy(6 + 4  * nu * ((-3) + 2  * nu), lr2 + x3 + y3)))
        end

        function J1312d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * ((
                        - 1) + 2  * nu) * (-x1 + y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) + lr1 ^ (-1) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1)
                         * (x3 - y3) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 -
            y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + lr2 ^ (-1) * (x1 - y1) * (lr2 +
            x2 - y2) ^ (-1) * ((7 + 8  * ((-2) + nu) * nu) * x3 + y3 + 8  * ((-1) + nu) *
            nu * y3) -4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 +
            y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * (
            (-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) ^ 3  * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) * (x2 - y2) *
            y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * atan((x1 - y1) * (x2 - y2) ^ (-1)) + 4  * ((-1) + nu) * ((
                        - 1) +2  * nu) * atan(lr2 * (x1 - y1), (x2 - y2) * (x3 + y3))))
        end

        function J1312d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * ((
                        - 1) + 2  * nu) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         ^ (-1) + (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + lr1 ^ (-1) * (x2 -
            y2) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 ^ (-1) * (x2 - y2) ^ 2  * (lr2 + x3 + y3) ^ (-1) + (lr2 + x2 - y2)
                         ^ (-1) * ((7 + 8  * ((-2) + nu) * nu) * x3 + y3 + 8  * ((-1) + nu) * nu * y3) + lr2 ^ (
            -1) * (x2 - y2) * (lr2 + x2 - y2) ^ (-1) * ((7 + 8  * ((-2) + nu) * nu)
                         * x3 + y3 + 8  * ((-1) + nu) * nu * y3) -4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (
            x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 +
            y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * x3 * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) - 2  * x3 * (x2 - y2) ^ 2  * y3 * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) + xlogy(4  * ((-1) + nu) * ((-1) + 2  * nu), lr2 + x3 + y3)))
        end

        function J1312d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x2 + (
            -1) * y2) ^ (-1) * (x3 - y3) ^ 2 + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x2 + (
            -1) * y2) * (lr2 + x3 + y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (
            x2 - y2) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + lr2 ^ (-1) * (lr2 + x2 - y2)
                         ^ (-1) * (x3 + y3) * ((7 + 8  * ((-2) + nu) * nu) * x3 + y3 + 8  * ((-1) + nu) * nu *
            y3) -4  * lr2 ^ (-1) * x3 * (x2 - y2) * y3 * (x3 + y3) ^ 2  * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 + (
            -1) * y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * (x2 - y2) * y3 * (2  * x3 + y3) * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x2 - y2) * y3 * (x3 + y3) ^ 2  *
            ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + log(lr1 + x2 - y2) + xlogy(7 + 8  * ((-2) +
                nu) * nu, lr2 + x2 - y2)))
        end

        function J1313d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x1 - y1) + 2  *
            (7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) - 2  * lr2 ^ (-1) * (x1 +
            (-1) * y1) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) -2  * ((-3) + 4  * nu) * x3 * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + (x1 - y1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (((-7) - 8  * ((-2) + nu) * nu) *
            x1 ^ 2 + 2  * (7 + 8  * ((-2) + nu) * nu) * x1 * y1 + ((-7) - 8  * ((-2) + nu) * nu) *
            y1 ^ 2 - 7  * (x3 ^ 2 + (x2 - y2) ^ 2) - 16  * x3 * y3 - 7  * y3 ^ 2 + (
            -8) * ((-2) + nu) * nu * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2))) * G ^ (-1))
        end

        function J1313d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x2 - y2) + 2  *
            (7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (x2 - y2) - 2  * lr2 ^ (-1) * (x2 +
            (-1) * y2) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) -2  * ((-3) + 4  * nu) * x3 * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + (x2 - y2) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (((-7) - 8  * ((-2) + nu) * nu) *
            x1 ^ 2 + 2  * (7 + 8  * ((-2) + nu) * nu) * x1 * y1 + ((-7) - 8  * ((-2) + nu) * nu) *
            y1 ^ 2 - 7  * (x3 ^ 2 + (x2 - y2) ^ 2) - 16  * x3 * y3 - 7  * y3 ^ 2 + (
            -8) * ((-2) + nu) * nu * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2))) * G ^ (-1))
        end

        function J1313d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x3 + (-1)
                         * y3) +2  * lr2 ^ (-1) * ((7 + 8  * ((-2) + nu) * nu) * x3 + 8  * ((-1) + nu) ^ 2  * y3) +
            2  * (lr2 + x3 + y3) ^ (-1) * ((-3) * x3 - 2  * y3 + 6  * nu * (x3 + y3) - 4  *
            nu ^ 2  * (x3 + y3)) -2  * lr2 ^ (-1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 +
            2  * y3 - 6  * nu * (x3 + y3) + 4  * nu ^ 2  * (x3 + y3)) -2  * ((-3) + 4  * nu) *
            x3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) ^ 2  * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + (x3 + y3) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (((-7) - 8  * ((
                        - 2) + nu) * nu) * x1 ^ 2 + 2  * (7 + 8  * ((-2) + nu) * nu) * x1 * y1 + ((-7) - 8  * ((
                        - 2) + nu) * nu) * y1 ^ 2 - 7  * (x3 ^ 2 + (x2 - y2) ^ 2) - 16  * x3 * y3 + (
            -7) * y3 ^ 2 - 8  * ((-2) + nu) * nu * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)) +2  *
            ((-3) + 4  * nu) * lr2 ^ (-1) * x3 * (1 - (x3 + y3) ^ 2  * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)) ^ (-1) + 2  * ((-3) + 4  * nu) *
            real(acoth(complex(lr2 ^ (-1) * (x3 + y3)))) + xlogy((-6) + 4  * (3 - 2  * nu) * nu, lr2 + x3 + y3)))
        end

        function J1323d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x1 + (-1)
                         * y1) ^ 2  * (lr1 + x2 - y2) ^ (-1) + (7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (
            x1 - y1) ^ 2  * (lr2 + x2 - y2) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu)
                         * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 +
            y3) +4  * lr2 ^ (-1) * x3 * (x1 - y1) ^ 2  * (x2 - y2) * y3 * ((x1 + (
            -1) * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * x3 * (-x2 + y2) * y3 *
            ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((-3) *
            x3 - y3 + 2  * nu * (x3 + y3)) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + (
            -4) * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 +
            (-1) * y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((-3) * x3 - y3 +
            2  * nu * (x3 + y3)) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * x3 * (x1 + (
            -1) * y1) ^ 2  * (x2 - y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + log(
                lr1 + x2 - y2) + xlogy(7 + 8  * ((-2) + nu) * nu, lr2 + x2 - y2)))
        end

        function J1323d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((x1 - y1) * (lr1 + x2 + (-1)
                         * y2) ^ (-1) + lr1 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr1 + x2 -
            y2) ^ (-1) + (7 + 8  * ((-2) + nu) * nu) * (x1 - y1) * (lr2 + x2 - y2) ^ (
            -1) +(7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) *
            (lr2 + x2 - y2) ^ (-1) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) + 2  * lr2 ^ (
            -1) * x3 * (-x1 + y1) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            4  * ((-1) + nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x3 + y3) * ((-3) * x3 - y3 + 2  * nu * (x3 + y3)) * ((x1 + (-1)
                         * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  *
            (x3 + y3) * ((-3) * x3 - y3 + 2  * nu * (x3 + y3)) * ((x1 - y1) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) + 2  * x3 * (x1 - y1) * (x2 - y2) ^ 2  * y3 * ((x1 +
            (-1) * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2)) * G ^ (-1))
        end

        function J1323d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x1 + (-1)
                         * y1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + (7 + 8  * ((-2) + nu) * nu)
                         * lr2 ^ (-1) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) + 4  *
            lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 + (-1)
                         * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) - 2  * lr2 ^ (-1) * (x1 - y1) * (x2 + (-1)
                         * y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) * ((-3) * x3 - y3 + 2  * nu * (x3 + y3)) * ((x1 - y1) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * ((-3)
                         * x3 - y3 + 2  * nu * (x3 + y3)) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1)
            +2  * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1) ^ (
                -1) * (x2 - y2)) -4  * ((-1) + nu) * ((-3) + 2  * nu) * atan(lr2 * (x1 + (
                    -1) * y1), (x2 - y2) * (x3 + y3))))
        end

        function J2312d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * ((
                        - 1) + 2  * nu) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) ^ 2 + (lr1 + x1 - y1) ^ (-1) * (x3 - y3) + lr1 ^ (-1) * (x1 -
            y1) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x3 + y3) ^ (-1) + (lr2 + x1 - y1)
                         ^ (-1) * ((7 + 8  * ((-2) + nu) * nu) * x3 + y3 + 8  * ((-1) + nu) * nu * y3) + lr2 ^ (
            -1) * (x1 - y1) * (lr2 + x1 - y1) ^ (-1) * ((7 + 8  * ((-2) + nu) * nu)
                         * x3 + y3 + 8  * ((-1) + nu) * nu * y3) -4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (
            (x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 +
            y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * x3 * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) ^ 2  * y3 * (x3 + y3) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) + xlogy(4  * ((-1) + nu) * ((-1) + 2  * nu), lr2 + x3 + y3)))
        end

        function J2312d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * ((-4) * ((-1) + nu) * (
            (-1) + 2  * nu) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x2 - y2) + lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2)
                         * (x3 - y3) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 -
            y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + lr2 ^ (-1) * (lr2 + x1 - y1) ^ (
            -1) * (x2 - y2) * ((7 + 8  * ((-2) + nu) * nu) * x3 + y3 + 8  * ((-1) + nu) *
            nu * y3) -4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 +
            y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x2 - y2) * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * (
            (-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 3  * (x3 + y3) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) * (x2 - y2)
                         * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) - 4  * ((-1) + nu) * ((-1) +
            2  * nu) * atan((x1 - y1) ^ (-1) * (x2 - y2)) + 4  * ((-1) + nu) * ((
                        - 1) +2  * nu) * atan(lr2 * (x2 - y2), (x1 - y1) * (x3 + y3))))
        end

        function J2312d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x3 - y3) ^ 2 + 4  * ((-1) + nu) * ((-1) + 2  * nu) * (x1 + (
            -1) * y1) * (lr2 + x3 + y3) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 ^ (-1) * (
            x1 - y1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) + lr2 ^ (-1) * (lr2 + x1 - y1)
                         ^ (-1) * (x3 + y3) * ((7 + 8  * ((-2) + nu) * nu) * x3 + y3 + 8  * ((-1) + nu) * nu *
            y3) -4  * lr2 ^ (-1) * x3 * (x1 - y1) * y3 * (x3 + y3) ^ 2  * ((x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * lr2 * (x1 + (
            -1) * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2)
                         ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * ((-1) + 2  *
            nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2)
                         ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * (x1 - y1) * y3 * (2  * x3 + y3) * ((x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) * y3 * (x3 + y3) ^ 2  *
            ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + log(lr1 + x1 - y1) + xlogy(7 + 8  * ((-2) +
                nu) * nu, lr2 + x1 - y1)))
        end

        function J2313d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * ((lr1 + x1 - y1) ^ (-1) * (
            x2 - y2) +lr1 ^ (-1) * (x1 - y1) * (lr1 + x1 - y1) ^ (-1) * (x2 +
            (-1) * y2) +(7 + 8  * ((-2) + nu) * nu) * (lr2 + x1 - y1) ^ (-1) * (x2 -
            y2) +(7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) * (lr2 + x1 -
            y1) ^ (-1) * (x2 - y2) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) + 2  * lr2 ^ (
            -1) * x3 * (-x2 + y2) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) +
            4  * ((-1) + nu) * lr2 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * (x3 + y3) * ((-3) * x3 - y3 + 2  * nu * (x3 + y3)) * ((x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1)
                         ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (
            x3 + y3) * ((-3) * x3 - y3 + 2  * nu * (x3 + y3)) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) + 2  * x3 * (x1 - y1) ^ 2  * (x2 - y2) * y3 * ((x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 +
            (x3 + y3) ^ 2) ^ (-3 / 2)) * G ^ (-1))
        end

        function J2313d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) ^ 2 + (7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (
            lr2 + x1 - y1) ^ (-1) * (x2 - y2) ^ 2 + 4  * ((-1) + nu) * ((-1) + 2  * nu)
                         * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 +
            y3) +4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) ^ 2  * y3 * ((x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) + 2  * lr2 ^ (-1) * x3 * (-x1 + y1) * y3 *
            ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 * (x1 + (-1)
                         * y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((-3) *
            x3 - y3 + 2  * nu * (x3 + y3)) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + (
            -4) * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((-3) * x3 - y3 +
            2  * nu * (x3 + y3)) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * x3 * (x1 + (
            -1) * y1) * (x2 - y2) ^ 2  * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + log(
                lr1 + x1 - y1) + xlogy(7 + 8  * ((-2) + nu) * nu, lr2 + x1 - y1)))
        end

        function J2313d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) * (x3 - y3) + (7 + 8  * ((-2) + nu) * nu)
                         * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 + y3) + 4  *
            lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 + (-1)
                         * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) - 2  * lr2 ^ (-1) * (x1 - y1) * (x2 + (-1)
                         * y2) * y3 * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 4  * ((-1) + nu) * lr2 * (
            x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) * ((-3) * x3 - y3 + 2  * nu * (x3 + y3)) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 4  * ((-1) + nu) * lr2 ^ (-1) * (x1 - y1) * ((x1 + (-1)
                         * y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * ((-3)
                         * x3 - y3 + 2  * nu * (x3 + y3)) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)
            +2  * x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 - y2)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) - 4  * ((-1) + nu) * ((-1) + 2  * nu) * atan((x1 - y1) * (
                x2 - y2) ^ (-1)) -4  * ((-1) + nu) * ((-3) + 2  * nu) * atan(lr2 * (x2 + (
                    -1) * y2), (x1 - y1) * (x3 + y3))))
        end
        function J2323d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x1 - y1) + 2  *
            (7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (x1 - y1) - 2  * lr2 ^ (-1) * (x1 +
            (-1) * y1) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) -2  * ((-3) + 4  * nu) * x3 * (x1 - y1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + (x1 - y1) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (((-7) - 8  * ((-2) + nu) * nu) *
            x1 ^ 2 + 2  * (7 + 8  * ((-2) + nu) * nu) * x1 * y1 + ((-7) - 8  * ((-2) + nu) * nu) *
            y1 ^ 2 - 7  * (x3 ^ 2 + (x2 - y2) ^ 2) - 16  * x3 * y3 - 7  * y3 ^ 2 + (
            -8) * ((-2) + nu) * nu * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2))) * G ^ (-1))
        end

        function J2323d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * (lr1 ^ (-1) * (x2 - y2) + 2  *
            (7 + 8  * ((-2) + nu) * nu) * lr2 ^ (-1) * (x2 - y2) - 2  * lr2 ^ (-1) * (x2 +
            (-1) * y2) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 + 2  * y3 - 6  * nu * (x3 + y3) + 4  *
            nu ^ 2  * (x3 + y3)) -2  * ((-3) + 4  * nu) * x3 * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + (x2 - y2) * ((x1 - y1) ^ 2 + (
            x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (((-7) - 8  * ((-2) + nu) * nu) *
            x1 ^ 2 + 2  * (7 + 8  * ((-2) + nu) * nu) * x1 * y1 + ((-7) - 8  * ((-2) + nu) * nu) *
            y1 ^ 2 - 7  * (x3 ^ 2 + (x2 - y2) ^ 2) - 16  * x3 * y3 - 7  * y3 ^ 2 + (
            -8) * ((-2) + nu) * nu * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2))) * G ^ (-1))
        end

        function J2323d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((-1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (lr1 ^ (-1) * (x3 + (-1)
                         * y3) +2  * lr2 ^ (-1) * ((7 + 8  * ((-2) + nu) * nu) * x3 + 8  * ((-1) + nu) ^ 2  * y3) +
            2  * (lr2 + x3 + y3) ^ (-1) * ((-3) * x3 - 2  * y3 + 6  * nu * (x3 + y3) - 4  *
            nu ^ 2  * (x3 + y3)) -2  * lr2 ^ (-1) * (x3 + y3) * (lr2 + x3 + y3) ^ (-1) * (3  * x3 +
            2  * y3 - 6  * nu * (x3 + y3) + 4  * nu ^ 2  * (x3 + y3)) -2  * ((-3) + 4  * nu) *
            x3 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x3 + y3) ^ 2  * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1 / 2) + (x3 + y3) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * (((-7) - 8  * ((
                        - 2) + nu) * nu) * x1 ^ 2 + 2  * (7 + 8  * ((-2) + nu) * nu) * x1 * y1 + ((-7) - 8  * ((
                        - 2) + nu) * nu) * y1 ^ 2 - 7  * (x3 ^ 2 + (x2 - y2) ^ 2) - 16  * x3 * y3 + (
            -7) * y3 ^ 2 - 8  * ((-2) + nu) * nu * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)) +2  *
            ((-3) + 4  * nu) * lr2 ^ (-1) * x3 * (1 - (x3 + y3) ^ 2  * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1)) ^ (-1) + 2  * ((-3) + 4  * nu) *
            real(acoth(complex(lr2 ^ (-1) * (x3 + y3)))) + xlogy((-6) + 4  * (3 - 2  * nu) * nu, lr2 + x3 + y3)))
        end

        function J3312d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (9  * x3 ^ 2  * (9  *
            x3 ^ 2 + (x1 - y1) ^ 2) ^ (-1) + 4  * nu ^ 2  * x3 ^ 2  * (nu ^ 2  * x3 ^ 2 + (x1 + (
            -1) * y1) ^ 2) ^ (-1) - ((-3) + 4  * nu) * (lr1 + x1 - y1) ^ (-1) * (
            x2 - y2) -((-3) + 4  * nu) * lr1 ^ (-1) * (x1 - y1) * (lr1 + x1 + (
            -1) * y1) ^ (-1) * (x2 - y2) + (5 + 4  * nu * ((-3) + 2  * nu)) * (lr2 + x1 + (-1)
                         * y1) ^ (-1) * (x2 - y2) + (5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 +
            (-1) * y1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) - ((-3) + 4  *
            nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * (lr1 + x2 - y2) ^ (-1) + (5 + 4  * nu *
            ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x2 - y2) ^ (-1) +
            2  * ((-1) + 2  * nu) * lr1 * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 -
            y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 + (-1)
                         * y3) ^ 2 - 2  * ((-1) + 2  * nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * (x2 + (-1)
                         * y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2)
                         ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 + 9  * y3 ^ 2  * ((x1 + (-1)
                         * y1) ^ 2 + 9  * y3 ^ 2) ^ (-1) + 4  * nu ^ 2  * y3 ^ 2  * ((x1 - y1) ^ 2 +
            nu ^ 2  * y3 ^ 2) ^ (-1) - 2  * (1 - 2  * nu) ^ 2  * lr2 * (x2 - y2) * (x3 +
            y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) ^ 2  * (x2 -
            y2) * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) * ((x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3)
                         ^ 2) -2  * x3 * (x1 - y1) ^ 2  * (x2 - y2) * y3 * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3) ^ 2) +2  * lr2 ^ (-1) * (x2 - y2) * ((
                x1- y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-1) * (x1 ^ 2  * x3 ^ 2 - 2  * x1 * x3 ^ 2  * y1 + x3 ^ 2  * y1 ^ 2 + 5  * x1 ^ 2  *
            x3 * y3 + x2 ^ 2  * x3 * y3 + 2  * x3 ^ 3  * y3 - 10  * x1 * x3 * y1 * y3 + 5  * x3 *
            y1 ^ 2  * y3 - 2  * x2 * x3 * y2 * y3 + x3 * y2 ^ 2  * y3 + x1 ^ 2  * y3 ^ 2 + 4  *
            x3 ^ 2  * y3 ^ 2 - 2  * x1 * y1 * y3 ^ 2 + y1 ^ 2  * y3 ^ 2 + 2  * x3 * y3 ^ 3 - 4  *
            nu * (x1 - y1) ^ 2  * (x3 + y3) ^ 2 + 4  * nu ^ 2  * (x1 - y1) ^ 2  * (x3 +
            y3) ^ 2) +xlogy(3 - 4  * nu, lr1 + x2 - y2) + xlogy(5 + 4  * nu * ((-3) + 2  *
                nu), lr2 + x2 - y2)))
        end

        function J3312d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (25  * x3 ^ 2  * (25  *
            x3 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 36  * nu ^ 2  * x3 ^ 2  * (9  * nu ^ 2  * x3 ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) + 8  * nu ^ 4  * x3 ^ 2  * (nu ^ 4  * x3 ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) - ((-3) + 4  * nu) * lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1)
                         * (x2 - y2) ^ 2 + (5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (lr2 + x1 + (-1)
                         * y1) ^ (-1) * (x2 - y2) ^ 2 - ((-3) + 4  * nu) * (x1 - y1) * (
            lr1 + x2 - y2) ^ (-1) - ((-3) + 4  * nu) * lr1 ^ (-1) * (x1 - y1)
                         * (x2 - y2) * (lr1 + x2 - y2) ^ (-1) + (5 + 4  * nu * ((-3) + 2  * nu)) * (
            x1 - y1) * (lr2 + x2 - y2) ^ (-1) + (5 + 4  * nu * ((-3) + 2  * nu)) *
            lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x2 - y2) ^ (-1) + 2  *
            ((-1) + 2  * nu) * lr1 * (x1 - y1) * ((x1 - y1) ^ 2 + (x3 - y3)
                         ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 -
            y3) ^ 2 - 2  * ((-1) + 2  * nu) * lr1 ^ (-1) * (x1 - y1) * (x2 - y2)
                         ^ 2  * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 +
            (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) ^ 2 + 25  * y3 ^ 2  * ((x2 -
            y2) ^ 2 + 25  * y3 ^ 2) ^ (-1) + 36  * nu ^ 2  * y3 ^ 2  * ((x2 - y2) ^ 2 + 9  *
            nu ^ 2  * y3 ^ 2) ^ (-1) + 8  * nu ^ 4  * y3 ^ 2  * ((x2 - y2) ^ 2 + nu ^ 4  *
            y3 ^ 2) ^ (-1) - 2  * (1 - 2  * nu) ^ 2  * lr2 * (x1 - y1) * (x3 + y3)
                         ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 4  * lr2 ^ (-1) * x3 * (x1 - y1) * (x2 - y2)
                         ^ 2  * y3 * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (
            x3 + y3) ^ 2) ^ (-2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3)
                         ^ 2) -2  * x3 * (x1 - y1) * (x2 - y2) ^ 2  * y3 * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3) ^ 2) +2  * lr2 ^ (-1) * (x1 - y1) * ((
                x1- y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2)
                         ^ (-1) * (x2 ^ 2  * x3 ^ 2 - 2  * x2 * x3 ^ 2  * y2 + x3 ^ 2  * y2 ^ 2 + x1 ^ 2  * x3 *
            y3 + 5  * x2 ^ 2  * x3 * y3 + 2  * x3 ^ 3  * y3 - 2  * x1 * x3 * y1 * y3 + x3 * y1 ^ 2  *
            y3 - 10  * x2 * x3 * y2 * y3 + 5  * x3 * y2 ^ 2  * y3 + x2 ^ 2  * y3 ^ 2 + 4  * x3 ^ 2  *
            y3 ^ 2 - 2  * x2 * y2 * y3 ^ 2 + y2 ^ 2  * y3 ^ 2 + 2  * x3 * y3 ^ 3 - 4  * nu * (x2 +
            (-1) * y2) ^ 2  * (x3 + y3) ^ 2 + 4  * nu ^ 2  * (x2 - y2) ^ 2  * (x3 + y3) ^ 2) +
            xlogy(3 - 4  * nu, lr1 + x1 - y1) + xlogy(5 + 4  * nu * ((-3) + 2  * nu), lr2 +
                x1 - y1)))
        end

        function J3312d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (9  * x3 * (9  * x3 ^ 2 + (
            x1 - y1) ^ 2) ^ (-1) * (-x1 + y1) + 4  * nu ^ 2  * x3 * (nu ^ 2  * x3 ^ 2 +
            (x1 - y1) ^ 2) ^ (-1) * (-x1 + y1) + 25  * x3 * (25  * x3 ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (-x2 + y2) + 36  * nu ^ 2  * x3 * (9  * nu ^ 2  * x3 ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (-x2 + y2) + 8  * nu ^ 4  * x3 * (nu ^ 4  * x3 ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * (-x2 + y2) - ((-3) + 4  * nu) * lr1 ^ (
            -1) * (lr1 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 - y3) - ((
                        - 3) +4  * nu) * lr1 ^ (-1) * (x1 - y1) * (lr1 + x2 - y2) ^ (-1) * (x3 + (
            -1) * y3) -2  * ((-1) + 2  * nu) * lr1 * (x1 - y1) * (x2 - y2) * ((
                x1- y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + (
            -1) * y3) ^ 2) ^ (-1) * (x3 - y3) - 2  * ((-1) + 2  * nu) * lr1 ^ (-1) * (
            x1 - y1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3)
                         ^ 3 + (5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (
            x2 - y2) * (x3 + y3) + (5 + 4  * nu * ((-3) + 2  * nu)) * lr2 ^ (-1) * (x1 + (-1)
                         * y1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) + 2  * (1 - 2  * nu) ^ 2  * lr2 * (
            x1 - y1) * (x2 - y2) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 4  * lr2 ^ (-1) *
            x3 * (x1 - y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (
            x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3) ^ 2) -4  * lr2 ^ (-1) * x3 * (x1 + (
            -1) * y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-2) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3) ^ 2) -2  * x3 * (x1 - y1) * (x2 + (
            -1) * y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 +
            (x3 + y3) ^ 2) ^ (-3 / 2) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + 2  * (x3 + y3)
                         ^ 2) +2  * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * ((x1 - y1) ^ 2 +
            (x3 + y3) ^ 2) ^ (-1) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * (x3 ^ 3 + (
            x1 ^ 2 + x2 ^ 2) * y3 + 9  * x3 ^ 2  * y3 - 2  * x1 * y1 * y3 + y1 ^ 2  * y3 - 2  *
            x2 * y2 * y3 + y2 ^ 2  * y3 + 11  * x3 * y3 ^ 2 + 3  * y3 ^ 3 - 4  * nu * (x3 + y3) ^ 3 +
            4  * nu ^ 2  * (x3 + y3) ^ 3) -3  * atan(3  * x3, x1 - y1) - 5  * atan(
                5  * x3, x2 - y2) + 12  * nu * atan((-3) * nu * x3, x2 - y2) + 4  * nu *
            atan(-nu * x3, x1 - y1) - 8  * nu ^ 2  * atan(nu ^ 2  * x3, x2 + (
                -1) * y2) +((-2) + 4  * nu) * atan(lr1 * (-x3 + y3), (x1 - y1) * (x2 +
                    (-1) * y2)) +2  * (1 - 2  * nu) ^ 2  * atan(lr2 * (x3 + y3), (x1 - y1) * (
                        x2 - y2))))
        end

        function J3313d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (9  * x2 ^ 2  * (9  *
            x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) + 25  * x2 ^ 2  * (25  * x2 ^ 2 + (x1 - y1)
                         ^ 2) ^ (-1) + 4  * nu ^ 2  * x2 ^ 2  * (nu ^ 2  * x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) +
            36  * nu ^ 2  * x2 ^ 2  * (9  * nu ^ 2  * x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) + 8  *
            nu ^ 4  * x2 ^ 2  * (nu ^ 4  * x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) + 4  * x3 * (lr2 + (
            -1) * x1 + y1) ^ (-1) + 4  * lr2 ^ (-1) * x3 * (-x1 + y1) * (lr2 - x1 + y1)
                         ^ (-1) + 9  * y2 ^ 2  * ((x1 - y1) ^ 2 + 9  * y2 ^ 2) ^ (-1) + 25  * y2 ^ 2  * ((
                x1- y1) ^ 2 + 25  * y2 ^ 2) ^ (-1) + 4  * nu ^ 2  * y2 ^ 2  * ((x1 - y1)
                         ^ 2 + nu ^ 2  * y2 ^ 2) ^ (-1) + 36  * nu ^ 2  * y2 ^ 2  * ((x1 - y1) ^ 2 + 9  *
            nu ^ 2  * y2 ^ 2) ^ (-1) + 8  * nu ^ 4  * y2 ^ 2  * ((x1 - y1) ^ 2 + nu ^ 4  *
            y2 ^ 2) ^ (-1) - ((-3) + 4  * nu) * (lr1 + x1 - y1) ^ (-1) * (x3 + (-1)
                         * y3) -((-3) + 4  * nu) * lr1 ^ (-1) * (x1 - y1) * (lr1 + x1 -
            y1) ^ (-1) * (x3 - y3) + 4  * ((-1) + nu) * lr1 * ((x1 - y1) ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 -
            y3) ^ 2) ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 -
            y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2)
                         ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3)
            -4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * (lr1 + x3 - y3) ^ (-1) + (
            -8) * ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x1 - y1) ^ 2  * (lr2 + x3 + y3) ^ (-1) + (
            lr2 + x1 - y1) ^ (-1) * ((-7) * x3 - 5  * y3 + 12  * nu * (x3 + y3) - 8  *
            nu ^ 2  * (x3 + y3)) -lr2 ^ (-1) * (x1 - y1) * (lr2 + x1 - y1) ^ (
            -1) * (7  * x3 + 5  * y3 - 12  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) + 8  * ((-1) +
            nu) ^ 2  * lr2 * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 -
            y2) ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 8  * ((-1) +
            nu) ^ 2  * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * x3 * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-1) - 2  * x3 * (x1 - y1) ^ 2  * y3 * (x3 + y3) * ((x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 +
            y3) ^ 2) ^ (-3 / 2) + xlogy(4 - 4  * nu, lr1 + x3 - y3) + xlogy((-8) * ((-1)
                +nu) ^ 2, lr2 + x3 + y3)))
        end

        function J3313d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (9  * x2 * (9  * x2 ^ 2 + (
            x1 - y1) ^ 2) ^ (-1) * (-x1 + y1) + 25  * x2 * (25  * x2 ^ 2 + (x1 + (-1)
                         * y1) ^ 2) ^ (-1) * (-x1 + y1) + 4  * nu ^ 2  * x2 * (nu ^ 2  * x2 ^ 2 + (x1 + (
            -1) * y1) ^ 2) ^ (-1) * (-x1 + y1) + 36  * nu ^ 2  * x2 * (9  * nu ^ 2  * x2 ^ 2 +
            (x1 - y1) ^ 2) ^ (-1) * (-x1 + y1) + 8  * nu ^ 4  * x2 * (nu ^ 4  *
            x2 ^ 2 + (x1 - y1) ^ 2) ^ (-1) * (-x1 + y1) + 4  * lr2 ^ (-1) * x3 * (lr2 +
            (-1) * x1 + y1) ^ (-1) * (-x2 + y2) - ((-3) + 4  * nu) * lr1 ^ (-1) * (
            lr1 + x1 - y1) ^ (-1) * (x2 - y2) * (x3 - y3) - 4  * ((-1) +
            nu) * lr1 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x2 - y2) * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 + (
            -1) * y3) -4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 3  * ((x2 - y2) ^ 2 + (
            x3 - y3) ^ 2) ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * lr1 ^ (-1) * (
            x1 - y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) - 8  * ((-1) +
            nu) ^ 2  * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + (
            -1) * lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x2 - y2) * (7  * x3 + 5  *
            y3 - 12  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) -4  * lr2 ^ (-1) * x3 * (x1 + (
            -1) * y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-2) - 8  * ((-1) + nu) ^ 2  * lr2 * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x2 -
            y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 8  * ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x1 -
            y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 3  *
            (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 -
            y1) * (x2 - y2) * y3 * (x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + 5  *
            atan((-5) * x2, x1 - y1) - 3  * atan(3  * x2, x1 - y1) + 4  * nu *
            atan(-nu * x2, x1 - y1) - 12  * nu * atan(3  * nu * x2, x1 + (-1)
                         * y1) + 8  * nu ^ 2  * atan(-nu ^ 2  * x2, x1 - y1) + (4 - 4  * nu) *
            atan(lr1 * (x2 - y2), (x1 - y1) * (x3 - y3)) - 8  * ((-1) +
            nu) ^ 2  * atan(lr2 * (x2 - y2), (x1 - y1) * (x3 + y3))))
        end

        function J3313d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * lr1 *
            (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (-1)
                         * y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2) ^ (-1) - ((-3) +
            4  * nu) * lr1 ^ (-1) * (lr1 + x1 - y1) ^ (-1) * (x3 - y3) ^ 2 - 4  * (
            (-1) + nu) * lr1 ^ (-1) * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) ^ 2 - 4  * ((-1) + nu) * (x1 - y1) * (lr1 + x3 + (
            -1) * y3) ^ (-1) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) * (x3 + (-1)
                         * y3) * (lr1 + x3 - y3) ^ (-1) - 4  * lr2 ^ (-1) * x3 * (lr2 - x1 + y1)
                         ^ (-1) * (x3 + y3) - 8  * ((-1) + nu) ^ 2  * (x1 - y1) * (lr2 + x3 + y3) ^ (
            -1) -8  * ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x1 - y1) * (x3 + y3) * (lr2 + x3 +
            y3) ^ (-1) - lr2 ^ (-1) * (lr2 + x1 - y1) ^ (-1) * (x3 + y3) * (7  * x3 +
            5  * y3 - 12  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) -4  * lr2 ^ (-1) * x3 * (
            x1 - y1) * y3 * (x3 + y3) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-2) +
            8  * ((-1) + nu) ^ 2  * lr2 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) -8  * ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) ^ 2  * ((x2 + (
            -1) * y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * (x1 - y1) * y3 * (2  *
            x3 + y3) * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 -
            y1) * y3 * (x3 + y3) ^ 2  * ((x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 + (
            -1) * y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + xlogy((-4), lr2 + (
                -1) * x1 + y1) +xlogy(3 - 4  * nu, lr1 + x1 - y1) + xlogy((-7) + 4  * (3
                    -2  * nu) * nu, lr2 + x1 - y1)))
        end

        function J3323d1(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (9  * x1 * (9  * x1 ^ 2 + (
            x2 - y2) ^ 2) ^ (-1) * (-x2 + y2) + 25  * x1 * (25  * x1 ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (-x2 + y2) + 4  * nu ^ 2  * x1 * (nu ^ 2  * x1 ^ 2 + (x2 + (
            -1) * y2) ^ 2) ^ (-1) * (-x2 + y2) + 36  * nu ^ 2  * x1 * (9  * nu ^ 2  * x1 ^ 2 +
            (x2 - y2) ^ 2) ^ (-1) * (-x2 + y2) + 8  * nu ^ 4  * x1 * (nu ^ 4  *
            x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (-x2 + y2) + 4  * lr2 ^ (-1) * x3 * ((
                        - 1) * x1 + y1) * (lr2 - x2 + y2) ^ (-1) - ((-3) + 4  * nu) * lr1 ^ (-1) *
            (x1 - y1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) - 4  * ((-1) +
            nu) * lr1 * (x1 - y1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1)
                         * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 + (
            -1) * y3) -4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 3  * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (
            x3 - y3) ^ 2) ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * lr1 ^ (-1) * (
            x1 - y1) * (x2 - y2) * (lr1 + x3 - y3) ^ (-1) - 8  * ((-1) +
            nu) ^ 2  * lr2 ^ (-1) * (x1 - y1) * (x2 - y2) * (lr2 + x3 + y3) ^ (-1) + (
            -1) * lr2 ^ (-1) * (x1 - y1) * (lr2 + x2 - y2) ^ (-1) * (7  * x3 + 5  *
            y3 - 12  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) -4  * lr2 ^ (-1) * x3 * (x1 + (
            -1) * y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-2) - 8  * ((-1) + nu) ^ 2  * lr2 * (x1 - y1) * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) * ((x1 -
            y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 8  * ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x1 -
            y1) ^ 3  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) *
            (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x1 -
            y1) * (x2 - y2) * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + 5  *
            atan((-5) * x1, x2 - y2) - 3  * atan(3  * x1, x2 - y2) + 4  * nu *
            atan(-nu * x1, x2 - y2) - 12  * nu * atan(3  * nu * x1, x2 + (-1)
                         * y2) + 8  * nu ^ 2  * atan(-nu ^ 2  * x1, x2 - y2) + (4 - 4  * nu) *
            atan(lr1 * (x1 - y1), (x2 - y2) * (x3 - y3)) - 8  * ((-1) +
            nu) ^ 2  * atan(lr2 * (x1 - y1), (x2 - y2) * (x3 + y3))))
        end

        function J3323d2(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (9  * x1 ^ 2  * (9  *
            x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 25  * x1 ^ 2  * (25  * x1 ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) + 4  * nu ^ 2  * x1 ^ 2  * (nu ^ 2  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) +
            36  * nu ^ 2  * x1 ^ 2  * (9  * nu ^ 2  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 8  *
            nu ^ 4  * x1 ^ 2  * (nu ^ 4  * x1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 9  * y1 ^ 2  * (9  *
            y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 25  * y1 ^ 2  * (25  * y1 ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) + 4  * nu ^ 2  * y1 ^ 2  * (nu ^ 2  * y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) +
            36  * nu ^ 2  * y1 ^ 2  * (9  * nu ^ 2  * y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 8  *
            nu ^ 4  * y1 ^ 2  * (nu ^ 4  * y1 ^ 2 + (x2 - y2) ^ 2) ^ (-1) + 4  * x3 * (lr2 + (
            -1) * x2 + y2) ^ (-1) + 4  * lr2 ^ (-1) * x3 * (-x2 + y2) * (lr2 - x2 + y2)
                         ^ (-1) - ((-3) + 4  * nu) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) + (
            -1) * ((-3) + 4  * nu) * lr1 ^ (-1) * (x2 - y2) * (lr1 + x2 - y2) ^ (-1)
                         * (x3 - y3) + 4  * ((-1) + nu) * lr1 * (x1 - y1) ^ 2  * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x1 - y1)
                         ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) ^ 2  *
            ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) * (x3 - y3) - 4  * ((
                        - 1) +nu) * lr1 ^ (-1) * (x2 - y2) ^ 2  * (lr1 + x3 - y3) ^ (-1) - 8  *
            ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x2 - y2) ^ 2  * (lr2 + x3 + y3) ^ (-1) + (lr2 + x2 +
            (-1) * y2) ^ (-1) * ((-7) * x3 - 5  * y3 + 12  * nu * (x3 + y3) - 8  * nu ^ 2  * (
            x3 + y3)) -lr2 ^ (-1) * (x2 - y2) * (lr2 + x2 - y2) ^ (-1) * (
            7  * x3 + 5  * y3 - 12  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) +8  * ((-1) + nu)
                         ^ 2  * lr2 * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (
            -1) * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 8  * ((-1) + nu)
                         ^ 2  * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2)
                         ^ 2) ^ (-1) * (x2 - y2) ^ 2  * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * x3 * y3 * (x3 + y3) * ((x1 - y1) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-1) - 2  * x3 * (x2 - y2) ^ 2  * y3 * (x3 + y3) * ((x1 - y1)
                         ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3)
                         ^ 2) ^ (-3 / 2) + xlogy(4 - 4  * nu, lr1 + x3 - y3) + xlogy((-8) * ((-1) +
                nu) ^ 2, lr2 + x3 + y3)))
        end

        function J3323d3(y1::R, y2::R, y3::R) where R
            lr1 = r1(y1, y2, y3)
            lr2 = r2(y1, y2, y3)
            ((1 / 16) * (1 - nu) ^ (-1) * 1 / π * G ^ (-1) * (4  * ((-1) + nu) * lr1 *
            (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 + (
            -1) * y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2) ^ (-1) - ((-3) +
            4  * nu) * lr1 ^ (-1) * (lr1 + x2 - y2) ^ (-1) * (x3 - y3) ^ 2 - 4  * (
            (-1) + nu) * lr1 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 -
            y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 - y3) ^ 2)
                         ^ (-1) * (x3 - y3) ^ 2 - 4  * ((-1) + nu) * (x2 - y2) * (lr1 + x3 + (
            -1) * y3) ^ (-1) - 4  * ((-1) + nu) * lr1 ^ (-1) * (x2 - y2) * (x3 + (-1)
                         * y3) * (lr1 + x3 - y3) ^ (-1) - 4  * lr2 ^ (-1) * x3 * (lr2 - x2 + y2)
                         ^ (-1) * (x3 + y3) - 8  * ((-1) + nu) ^ 2  * (x2 - y2) * (lr2 + x3 + y3) ^ (
            -1) -8  * ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x2 - y2) * (x3 + y3) * (lr2 + x3 +
            y3) ^ (-1) - lr2 ^ (-1) * (lr2 + x2 - y2) ^ (-1) * (x3 + y3) * (7  * x3 +
            5  * y3 - 12  * nu * (x3 + y3) + 8  * nu ^ 2  * (x3 + y3)) -4  * lr2 ^ (-1) * x3 * (
            x2 - y2) * y3 * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-2) +
            8  * ((-1) + nu) ^ 2  * lr2 * (x1 - y1) ^ 2  * ((x1 - y1) ^ 2 + (x2 + (-1)
                         * y2) ^ 2) ^ (-1) * (x2 - y2) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (
            -1) -8  * ((-1) + nu) ^ 2  * lr2 ^ (-1) * (x1 - y1) ^ 2  * ((x1 - y1)
                         ^ 2 + (x2 - y2) ^ 2) ^ (-1) * (x2 - y2) * (x3 + y3) ^ 2  * ((x1 + (-1)
                         * y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) + 2  * lr2 ^ (-1) * (x2 - y2) * y3 * (2  * x3 +
            y3) * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) - 2  * x3 * (x2 - y2) *
            y3 * (x3 + y3) ^ 2  * ((x1 - y1) ^ 2 + (x3 + y3) ^ 2) ^ (-1) * ((x1 -
            y1) ^ 2 + (x2 - y2) ^ 2 + (x3 + y3) ^ 2) ^ (-3 / 2) + xlogy((-4), lr2 -
                x2 + y2) + xlogy(3 - 4  * nu, lr1 + x2 - y2) + xlogy((-7) + 4  * (3 - 2  * nu)
                         * nu, lr2 + x2 - y2)))
        end

        function IU1d1(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J1123d1(y1, y2, y3)
            +2 * G * epsv12p * (J1223d1(y1, y2, y3) + J1113d1(y1, y2, y3))
            +2 * G * epsv13p * (J1323d1(y1, y2, y3) + J1112d1(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J1213d1(y1, y2, y3)
            +2 * G * epsv23p * (J1212d1(y1, y2, y3) + J1313d1(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J1312d1(y1, y2, y3))
        end

        function IU1d2(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J1123d2(y1, y2, y3)
            +2 * G * epsv12p * (J1223d2(y1, y2, y3) + J1113d2(y1, y2, y3))
            +2 * G * epsv13p * (J1323d2(y1, y2, y3) + J1112d2(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J1213d2(y1, y2, y3)
            +2 * G * epsv23p * (J1212d2(y1, y2, y3) + J1313d2(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J1312d2(y1, y2, y3))
        end

        function IU1d3(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J1123d3(y1, y2, y3)
            +2 * G * epsv12p * (J1223d3(y1, y2, y3) + J1113d3(y1, y2, y3))
            +2 * G * epsv13p * (J1323d3(y1, y2, y3) + J1112d3(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J1213d3(y1, y2, y3)
            +2 * G * epsv23p * (J1212d3(y1, y2, y3) + J1313d3(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J1312d3(y1, y2, y3))
        end

        function IU2d1(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J2123d1(y1, y2, y3)
            +2 * G * epsv12p * (J2223d1(y1, y2, y3) + J2113d1(y1, y2, y3))
            +2 * G * epsv13p * (J2323d1(y1, y2, y3) + J2112d1(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J2213d1(y1, y2, y3)
            +2 * G * epsv23p * (J2212d1(y1, y2, y3) + J2313d1(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J2312d1(y1, y2, y3))
        end

        function IU2d2(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J2123d2(y1, y2, y3)
            +2 * G * epsv12p * (J2223d2(y1, y2, y3) + J2113d2(y1, y2, y3))
            +2 * G * epsv13p * (J2323d2(y1, y2, y3) + J2112d2(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J2213d2(y1, y2, y3)
            +2 * G * epsv23p * (J2212d2(y1, y2, y3) + J2313d2(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J2312d2(y1, y2, y3))
        end

        function IU2d3(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J2123d3(y1, y2, y3)
            +2 * G * epsv12p * (J2223d3(y1, y2, y3) + J2113d3(y1, y2, y3))
            +2 * G * epsv13p * (J2323d3(y1, y2, y3) + J2112d3(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J2213d3(y1, y2, y3)
            +2 * G * epsv23p * (J2212d3(y1, y2, y3) + J2313d3(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J2312d3(y1, y2, y3))
        end

        function IU3d1(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J3123d1(y1, y2, y3)
            +2 * G * epsv12p * (J3223d1(y1, y2, y3) + J3113d1(y1, y2, y3))
            +2 * G * epsv13p * (J3323d1(y1, y2, y3) + J3112d1(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J3213d1(y1, y2, y3)
            +2 * G * epsv23p * (J3212d1(y1, y2, y3) + J3313d1(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J3312d1(y1, y2, y3))
        end

        function IU3d2(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J3123d2(y1, y2, y3)
            +2 * G * epsv12p * (J3223d2(y1, y2, y3) + J3113d2(y1, y2, y3))
            +2 * G * epsv13p * (J3323d2(y1, y2, y3) + J3112d2(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J3213d2(y1, y2, y3)
            +2 * G * epsv23p * (J3212d2(y1, y2, y3) + J3313d2(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J3312d2(y1, y2, y3))
        end

        function IU3d3(y1::R, y2::R, y3::R) where R
            ((lambda * epsvkk + 2 * G * epsv11p) * J3123d3(y1, y2, y3)
            +2 * G * epsv12p * (J3223d3(y1, y2, y3) + J3113d3(y1, y2, y3))
            +2 * G * epsv13p * (J3323d3(y1, y2, y3) + J3112d3(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv22p) * J3213d3(y1, y2, y3)
            +2 * G * epsv23p * (J3212d3(y1, y2, y3) + J3313d3(y1, y2, y3))
            +(lambda * epsvkk + 2 * G * epsv33p) * J3312d3(y1, y2, y3))
        end

        e11p =
            (IU1d1(L, T / 2, q3 + W) - IU1d1(L, -T / 2, q3 + W) + IU1d1(L, -T / 2, q3) - IU1d1(L, T / 2, q3)
            -IU1d1(zero(R), T / 2, q3 + W) + IU1d1(zero(R), -T / 2, q3 + W) - IU1d1(zero(R), -T / 2, q3) + IU1d1(zero(R), T / 2, q3))
        u12 =
            (IU1d2(L, T / 2, q3 + W) - IU1d2(L, -T / 2, q3 + W) + IU1d2(L, -T / 2, q3) - IU1d2(L, T / 2, q3)
            -IU1d2(zero(R), T / 2, q3 + W) + IU1d2(zero(R), -T / 2, q3 + W) - IU1d2(zero(R), -T / 2, q3) + IU1d2(zero(R), T / 2, q3))
        u13 =
            (IU1d3(L, T / 2, q3 + W) - IU1d3(L, -T / 2, q3 + W) + IU1d3(L, -T / 2, q3) - IU1d3(L, T / 2, q3)
            -IU1d3(zero(R), T / 2, q3 + W) + IU1d3(zero(R), -T / 2, q3 + W) - IU1d3(zero(R), -T / 2, q3) + IU1d3(zero(R), T / 2, q3))
        u21 =
            (IU2d1(L, T / 2, q3 + W) - IU2d1(L, -T / 2, q3 + W) + IU2d1(L, -T / 2, q3) - IU2d1(L, T / 2, q3)
            -IU2d1(zero(R), T / 2, q3 + W) + IU2d1(zero(R), -T / 2, q3 + W) - IU2d1(zero(R), -T / 2, q3) + IU2d1(zero(R), T / 2, q3))
        e22p =
            (IU2d2(L, T / 2, q3 + W) - IU2d2(L, -T / 2, q3 + W) + IU2d2(L, -T / 2, q3) - IU2d2(L, T / 2, q3)
            -IU2d2(zero(R), T / 2, q3 + W) + IU2d2(zero(R), -T / 2, q3 + W) - IU2d2(zero(R), -T / 2, q3) + IU2d2(zero(R), T / 2, q3))
        u23 =
            (IU2d3(L, T / 2, q3 + W) - IU2d3(L, -T / 2, q3 + W) + IU2d3(L, -T / 2, q3) - IU2d3(L, T / 2, q3)
            -IU2d3(zero(R), T / 2, q3 + W) + IU2d3(zero(R), -T / 2, q3 + W) - IU2d3(zero(R), -T / 2, q3) + IU2d3(zero(R), T / 2, q3))
        u31 =
            (IU3d1(L, T / 2, q3 + W) - IU3d1(L, -T / 2, q3 + W) + IU3d1(L, -T / 2, q3) - IU3d1(L, T / 2, q3)
            -IU3d1(zero(R), T / 2, q3 + W) + IU3d1(zero(R), -T / 2, q3 + W) - IU3d1(zero(R), -T / 2, q3) + IU3d1(zero(R), T / 2, q3))
        u32 =
            (IU3d2(L, T / 2, q3 + W) - IU3d2(L, -T / 2, q3 + W) + IU3d2(L, -T / 2, q3) - IU3d2(L, T / 2, q3)
            -IU3d2(zero(R), T / 2, q3 + W) + IU3d2(zero(R), -T / 2, q3 + W) - IU3d2(zero(R), -T / 2, q3) + IU3d2(zero(R), T / 2, q3))
        e33 =
            (IU3d3(L, T / 2, q3 + W) - IU3d3(L, -T / 2, q3 + W) + IU3d3(L, -T / 2, q3) - IU3d3(L, T / 2, q3)
            -IU3d3(zero(R), T / 2, q3 + W) + IU3d3(zero(R), -T / 2, q3 + W) - IU3d3(zero(R), -T / 2, q3) + IU3d3(zero(R), T / 2, q3))

        e12p = (u12 + u21) / 2
        e13p = (u13 + u31) / 2
        e23p = (u23 + u32) / 2

        e11 = (cosd(theta) * e11p - sind(theta) * e12p) * cosd(theta) - (cosd(theta) * e12p - sind(theta) * e22p) * sind(theta)
        e12 = (cosd(theta) * e11p - sind(theta) * e12p) * sind(theta) + (cosd(theta) * e12p - sind(theta) * e22p) * cosd(theta)
        e13 = cosd(theta) * e13p - sind(theta) * e23p
        e22 = (sind(theta) * e11p + cosd(theta) * e12p) * sind(theta) + (sind(theta) * e12p + cosd(theta) * e22p) * cosd(theta)
        e23 = sind(theta) * e13p + cosd(theta) * e23p

        epsv11 = (cosd(theta) * epsv11p - sind(theta) * epsv12p) * cosd(theta) - (cosd(theta) * epsv12p - sind(theta) * epsv22p) * sind(theta)
        epsv12 = (cosd(theta) * epsv11p - sind(theta) * epsv12p) * sind(theta) + (cosd(theta) * epsv12p - sind(theta) * epsv22p) * cosd(theta)
        epsv13 = cosd(theta) * epsv13p - sind(theta) * epsv23p
        epsv22 = (sind(theta) * epsv11p + cosd(theta) * epsv12p) * sind(theta) + (sind(theta) * epsv12p + cosd(theta) * epsv22p) * cosd(theta)
        epsv23 = sind(theta) * epsv13p + cosd(theta) * epsv23p
        epsv33 = epsv33p

        ϵ[1] = e11 - epsv11 * S(x1 / L) * omega(x2 / T) * S((x3 - q3) / W)
        ϵ[2] = e12 - epsv12 * S(x1 / L) * omega(x2 / T) * S((x3 - q3) / W)
        ϵ[3] = e13 - epsv13 * S(x1 / L) * omega(x2 / T) * S((x3 - q3) / W)
        ϵ[4] = e22 - epsv22 * S(x1 / L) * omega(x2 / T) * S((x3 - q3) / W)
        ϵ[5] = e23 - epsv23 * S(x1 / L) * omega(x2 / T) * S((x3 - q3) / W)
        ϵ[6] = e33 - epsv33 * S(x1 / L) * omega(x2 / T) * S((x3 - q3) / W)

        return nothing
    end

end

heaviside(x::T) where T = x ≤ zero(T) ? zero(T) : one(T)

xlogy(x::T, y::T) where T = isapprox(x, zero(T)) ? zero(T) : x * log(y)

xlogy(x, y) = xlogy(promote(x, y)...)

omega(x::T) where T = heaviside(x + 1/2) - heaviside(x - 1/2)

S(x::T) where T = omega(x - 1/2)
