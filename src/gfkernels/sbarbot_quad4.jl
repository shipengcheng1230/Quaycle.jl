using StatsFuns

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
function sbarbot_disp_quad4!(
    x1::R, x2::R, x3::R, q1::R, q2::R, q3::R,
    L::R, T::R, W::R, theta::R,
    epsv11p::R, epsv12p::R, epsv13p::R, epsv22p::R, epsv23p::R, epsv33p::R,
    G::R, nu::R,
    ) where R

    function J1112(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
        -1) -4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)^ (
          -1)* (x2 - y2)) -x3* atan2(x3, x1 - y1) - 3 * x3*
        atan2(3 * x3, x1 - y1) + 4 * nu* x3* atan2(-nu* x3, x1 -
          y1) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan2(r2(y1, y2, y3)* (-x1 + y1), (
            x2-y2)* (x3 + y3)) -4 * ((-1) + nu)* (x3 - y3)* atan2(r1(y1, y2, y3)* (
              x3 - y3), (x1 - y1)* (x2 - y2)) +3 * y3* atan2((-3)* y3,
                x1 - y1) - y3* atan2(y3, x1 - y1) - 4 * nu* y3* atan2(
                  nu* y3, x1 - y1) - 4 * ((-1) + nu)* (x3 + y3)* atan2(r2(y1, y2, y3)* (x3 + y3), (
                    x1-y1)* (x2 - y2)) +xlogy(-((-3) + 4 * nu)* (x1 -
                      y1), r1(y1, y2, y3) + x2 - y2) +xlogy((5 + 4 * nu* ((-3) + 2 * nu))* (x1 - y1),
                        r2(y1, y2, y3) + x2 - y2) + xlogy((-4)* ((-1) + nu)* (x2 - y2), r1(y1, y2, y3) + x1 -
                          y1) + xlogy((-4)* ((-1) + nu)* (x2 - y2), r2(y1, y2, y3) + x1 - y1))
    end

    function J1113(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* (x1 + (
        -1)* y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (-((-1) +
        nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
        y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
        ) +x2* atan2(-x2, x1 - y1) - 3 * x2* atan2(3 * x2, x1 -
          y1) + 4 * nu* x2* atan2(-nu* x2, x1 - y1) - 4 * ((-1) + nu)* (
        x2 - y2)* atan2(r1(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 - y3)
        ) + 4 * ((-1) + nu)* (x2 - y2)* atan2(r2(y1, y2, y3)* (x2 - y2), (x1 -
          y1)* (x3 + y3)) +3 * y2* atan2((-3)* y2, x1 - y1) - y2* atan2(
            y2, x1 - y1) - 4 * nu* y2* atan2(nu* y2, x1 - y1) + xlogy((-1)
        * ((-3) + 4 * nu)* (x1 - y1), r1(y1, y2, y3) + x3 - y3) + xlogy(-(3
              -6 * nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +xlogy((-4)* ((-1) + nu)* (
                x3 - y3), r1(y1, y2, y3) + x1 - y1) +xlogy(4 * ((-1) + nu)* (x3 + y3), r2(y1, y2, y3) + x1 + (
                  -1)* y1))
    end

    function J1123(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)* ((
          x1- y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* ((x1 + (-1)
        * y1)^ 2 + (x3 + y3)^ 2)^ (-1)* (x3* ((x3^ 2 + (x1 - y1)^ 2)* (
        x3^ 2 + (x1 - y1)^ 2 + (x2 - y2)^ 2) +x3* (3 * x3^ 2 + 2 * (x1 + (-1)
        * y1)^ 2 + (x2 - y2)^ 2)* y3 + 3 * x3^ 2 * y3^ 2 + x3* y3^ 3) -((
        - 1) +nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)
        ^ 2) +((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)* y3* (2 * x3 + y3)* ((x1 - y1)
        ^ 2 + (x3 + y3)^ 2)) +2 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)* atan((
          x1-y1)* (x2 - y2)^ (-1)) +x1* atan2(-x1, x2 - y2)
        -3 * x1* atan2(3 * x1, x2 - y2) + 4 * nu* x1* atan2(-nu* x1,
          x2 - y2) + 3 * y1* atan2((-3)* y1, x2 - y2) - y1* atan2(
            y1, x2 - y2) - 4 * nu* y1* atan2(nu* y1, x2 - y2) + 2 * ((-1) +
        2 * nu)* (x1 - y1)* atan2(r1(y1, y2, y3)* (-x1 + y1), (x2 - y2)* (x3 +
          (-1)* y3)) +2 * (1 - 2 * nu)^ 2 * (x1 - y1)* atan2(r2(y1, y2, y3)* (-x1 +
            y1), (x2 - y2)* (x3 + y3)) +xlogy((-2)* x3, r2(y1, y2, y3) - x2 + y2) + xlogy((
        -1)* ((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x3 - y3) +xlogy(-(3 + (
              -6)* nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +xlogy(-((-3) + 4 *
                nu)* (x3 - y3), r1(y1, y2, y3) + x2 - y2) +xlogy(-(5 + 4 * nu* ((-3) + 2 *
                  nu))* (x3 + y3), r2(y1, y2, y3) + x2 - y2))
    end

    function J2112(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (-r1(y1, y2, y3) + (1 + 8 * ((
        - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + xlogy((-4)* ((-1) + nu)* ((
        - 1) + 2 * nu)* (x3 + y3), r2(y1, y2, y3) + x3 + y3))
    end

    function J2113(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* ((x1 +
        (-1)* y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* (-((-1) +
        nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
        y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
        ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +
        xlogy(-x2 + y2, r1(y1, y2, y3) + x3 - y3))
    end

    function J2123(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* (x1 + (
        -1)* y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (-((-1) +
        nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
        y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
        ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +
        xlogy(-x1 + y1, r1(y1, y2, y3) + x3 - y3))
    end

    function J3112(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
        x3* (x2 - y2)* y3* (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
        -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)* atan((x1 - y1)
        * (x2 - y2)^ (-1)) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)*
        atan2(r2(y1, y2, y3)* (-x1 + y1), (x2 - y2)* (x3 + y3)) + xlogy((-4)* ((-1) +
          nu)* ((-1) + 2 * nu)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +xlogy(x3 - y3, r1(y1, y2, y3) +
            x2 - y2) + xlogy(-x3 - 7 * y3 - 8 * nu^ 2 * (x3 + y3) + 8 * nu* (
              x3 + 2 * y3), r2(y1, y2, y3) + x2 - y2))
    end

    function J3113(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (r1(y1, y2, y3) + ((-1) - 8 * ((
        - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + 2 * ((-3) + 4 * nu)* x3*
        acoth(r2(y1, y2, y3)^ (-1)* (x3 + y3)) + xlogy(2 * (3 * x3 + 2 * y3 - 6 * nu* (x3 + y3) +
          4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3))
    end

    function J3123(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
        -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)^ (-1)
        * (x2 - y2)) + 4 * ((-1) + 2 * nu)* (nu* x3 + ((-1) + nu)* y3)* atan2(
          r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) + xlogy(x1 - y1, r1(y1, y2, y3) + x2 +
            (-1)* y2) + xlogy(-(1 + 8 * ((-1) + nu)* nu)* (x1 - y1), r2(y1, y2, y3) + x2 + (
              -1)* y2))
    end

    function J1212(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (-r1(y1, y2, y3) + (1 + 8 * ((
        - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + xlogy((-4)* ((-1) + nu)* ((
        - 1) + 2 * nu)* (x3 + y3), r2(y1, y2, y3) + x3 + y3))
    end

    function J1213(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* ((x1 +
        (-1)* y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* (-((-1) +
        nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
        y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
        ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +
        xlogy(-x2 + y2, r1(y1, y2, y3) + x3 - y3))
    end

    function J1223(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* (x1 + (
        -1)* y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (-((-1) +
        nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
        y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
        ) +xlogy(-((-1) - 2 * nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +
        xlogy(-x1 + y1, r1(y1, y2, y3) + x3 - y3))
    end

    function J2212(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x1 - y1)* (x2 - y2)* y3* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (
        -1) -4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)* (
          x2 - y2)^ (-1)) -x3* atan2(x3, x1 - y1) - 3 * x3*
        atan2(3 * x3, x1 - y1) + 4 * nu* x3* atan2(-nu* x3, x1 -
          y1) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan2(r2(y1, y2, y3)* (-x2 + y2), (
            x1-y1)* (x3 + y3)) -4 * ((-1) + nu)* (x3 - y3)* atan2(r1(y1, y2, y3)* (
              x3 - y3), (x1 - y1)* (x2 - y2)) +3 * y3* atan2((-3)* y3,
                x1 - y1) - y3* atan2(y3, x1 - y1) - 4 * nu* y3* atan2(
                  nu* y3, x1 - y1) - 4 * ((-1) + nu)* (x3 + y3)* atan2(r2(y1, y2, y3)* (x3 + y3), (
                    x1-y1)* (x2 - y2)) +xlogy((-4)* ((-1) + nu)* (x1 - y1),
                      r1(y1, y2, y3) + x2 - y2) + xlogy((-4)* ((-1) + nu)* (x1 - y1), r2(y1, y2, y3) + x2 -
                        y2) + xlogy(-((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x1 - y1) + xlogy(
                          (5 + 4 * nu* ((-3) + 2 * nu))* (x2 - y2), r2(y1, y2, y3) + x1 - y1))
    end

    function J2213(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)* (
        x1 - y1)* ((x1 - y1)^ 2 + (x2 - y2)^ 2)^ (-1)* ((x2 + (-1)
        * y2)^ 2 + (x3 + y3)^ 2)^ (-1)* (x3* ((x3^ 2 + (x2 - y2)^ 2)* (
        x3^ 2 + (x1 - y1)^ 2 + (x2 - y2)^ 2) +x3* (3 * x3^ 2 + (x1 -
        y1)^ 2 + 2 * (x2 - y2)^ 2)* y3 + 3 * x3^ 2 * y3^ 2 + x3* y3^ 3) -(
        (-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)
        ^ 2) +((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)* y3* (2 * x3 + y3)* ((x2 - y2)
        ^ 2 + (x3 + y3)^ 2)) +2 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)* atan((
          x1-y1)^ (-1)* (x2 - y2)) +x2* atan2(-x2, x1 - y1)
        -3 * x2* atan2(3 * x2, x1 - y1) + 4 * nu* x2* atan2(-nu* x2,
          x1 - y1) + 3 * y2* atan2((-3)* y2, x1 - y1) - y2* atan2(
            y2, x1 - y1) - 4 * nu* y2* atan2(nu* y2, x1 - y1) + 2 * ((-1) +
        2 * nu)* (x2 - y2)* atan2(r1(y1, y2, y3)* (-x2 + y2), (x1 - y1)* (x3 +
          (-1)* y3)) +2 * (1 - 2 * nu)^ 2 * (x2 - y2)* atan2(r2(y1, y2, y3)* (-x2 +
            y2), (x1 - y1)* (x3 + y3)) +xlogy((-2)* x3, r2(y1, y2, y3) - x1 + y1) + xlogy((
        -1)* ((-3) + 4 * nu)* (x1 - y1), r1(y1, y2, y3) + x3 - y3) +xlogy(-(3 + (
              -6)* nu + 4 * nu^ 2)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +xlogy(-((-3) + 4 *
                nu)* (x3 - y3), r1(y1, y2, y3) + x1 - y1) +xlogy(-(5 + 4 * nu* ((-3) + 2 *
                  nu))* (x3 + y3), r2(y1, y2, y3) + x1 - y1))
    end

    function J2223(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* ((x1 +
        (-1)* y1)^ 2 + (x2 - y2)^ 2)^ (-1)* (x2 - y2)* (-((-1) +
        nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)^ 2 * (x3 + y3) + ((-1) + nu)* ((-1) + 2 * nu)* r2(y1, y2, y3)*
        y3* (2 * x3 + y3) + x3* ((x1 - y1)^ 2 + (x2 - y2)^ 2 + x3* (x3 + y3))
        ) +x1* atan2(-x1, x2 - y2) - 3 * x1* atan2(3 * x1, x2 -
          y2) + 4 * nu* x1* atan2(-nu* x1, x2 - y2) - 4 * ((-1) + nu)* (
        x1 - y1)* atan2(r1(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 - y3)
        ) + 4 * ((-1) + nu)* (x1 - y1)* atan2(r2(y1, y2, y3)* (x1 - y1), (x2 -
          y2)* (x3 + y3)) +3 * y1* atan2((-3)* y1, x2 - y2) - y1* atan2(
            y1, x2 - y2) - 4 * nu* y1* atan2(nu* y1, x2 - y2) + xlogy((-1)
        * ((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x3 - y3) + xlogy(-(3
              -6 * nu + 4 * nu^ 2)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) +xlogy((-4)* ((-1) + nu)* (
                x3 - y3), r1(y1, y2, y3) + x2 - y2) +xlogy(4 * ((-1) + nu)* (x3 + y3), r2(y1, y2, y3) + x2 + (
                  -1)* y2))
    end

    function J3212(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
        x3* (x1 - y1)* y3* (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (
        -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)* atan((x1 - y1)
        ^ (-1)* (x2 - y2)) + 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)*
        atan2(r2(y1, y2, y3)* (-x2 + y2), (x1 - y1)* (x3 + y3)) + xlogy((-4)* ((-1) +
          nu)* ((-1) + 2 * nu)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) +xlogy(x3 - y3, r1(y1, y2, y3) +
            x1 - y1) + xlogy(-x3 - 7 * y3 - 8 * nu^ 2 * (x3 + y3) + 8 * nu* (
              x3 + 2 * y3), r2(y1, y2, y3) + x1 - y1))
    end

    function J3213(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x1 - y1)* (x2 - y2)* y3* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (
        -1) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 - y1)* (x2 + (
          -1)* y2)^ (-1)) +4 * ((-1) + 2 * nu)* (nu* x3 + ((-1) + nu)* y3)* atan2(
            r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) + xlogy(x2 - y2, r1(y1, y2, y3) + x1 +
              (-1)* y1) + xlogy(-(1 + 8 * ((-1) + nu)* nu)* (x2 - y2), r2(y1, y2, y3) + x1 + (
                -1)* y1))
    end

    function J3223(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (r1(y1, y2, y3) + ((-1) - 8 * ((
        - 1) + nu)* nu)* r2(y1, y2, y3) - 2 * r2(y1, y2, y3)^ (-1)* x3* y3 + 2 * ((-3) + 4 * nu)* x3*
        acoth(r2(y1, y2, y3)^ (-1)* (x3 + y3)) + xlogy(2 * (3 * x3 + 2 * y3 - 6 * nu* (x3 + y3) +
          4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3))
    end

    function J1312(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x2 - y2)* y3* (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (-1) + (
        -4)* ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)* atan((x1 - y1)* (
          x2 - y2)^ (-1)) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x1 - y1)*
        atan2(r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) + xlogy(4 * ((-1) + nu)
        * ((-1) + 2 * nu)* (x2 - y2), r2(y1, y2, y3) + x3 + y3) + xlogy(x3 - y3, r1(y1, y2, y3) + x2 + (
          -1)* y2) +xlogy((7 + 8 * ((-2) + nu)* nu)* x3 + y3 + 8 * ((-1) + nu)* nu* y3,
            r2(y1, y2, y3) + x2 - y2))
    end

    function J1313(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (r1(y1, y2, y3) + r2(y1, y2, y3)^ (-1)* ((7 +
        8 * ((-2) + nu)* nu)* r2(y1, y2, y3)^ 2 + 2 * x3* y3) +2 * ((-3) + 4 * nu)* x3* acoth(
          r2(y1, y2, y3)^ (-1)* (x3 + y3)) + xlogy(2 * ((-3)* x3 - 2 * y3 + 6 * nu* (x3 + y3)
            -4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3))
    end

    function J1323(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
        x3* (x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)
        ^ 2)^ (-1) - 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 -
          y1)^ (-1)* (x2 - y2)) -4 * ((-1) + nu)* ((-3)* x3 - y3 + 2 *
        nu* (x3 + y3))* atan2(r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) +
        xlogy(x1 - y1, r1(y1, y2, y3) + x2 - y2) + xlogy((7 + 8 * ((-2) + nu)* nu)* (x1 +
          (-1)* y1), r2(y1, y2, y3) + x2 - y2))
    end

    function J2312(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x1 - y1)* y3* (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (-1) + (
        -4)* ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)* atan((x1 - y1)^ (
          -1)* (x2 - y2)) +4 * ((-1) + nu)* ((-1) + 2 * nu)* (x2 - y2)*
        atan2(r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) + xlogy(4 * ((-1) + nu)
        * ((-1) + 2 * nu)* (x1 - y1), r2(y1, y2, y3) + x3 + y3) + xlogy(x3 - y3, r1(y1, y2, y3) + x1 + (
          -1)* y1) +xlogy((7 + 8 * ((-2) + nu)* nu)* x3 + y3 + 8 * ((-1) + nu)* nu* y3,
            r2(y1, y2, y3) + x1 - y1))
    end

    function J2313(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* ((-2)* r2(y1, y2, y3)^ (-1)*
        x3* (x1 - y1)* (x2 - y2)* y3* ((x2 - y2)^ 2 + (x3 + y3)
        ^ 2)^ (-1) - 4 * ((-1) + nu)* ((-1) + 2 * nu)* (x3 + y3)* atan((x1 -
          y1)* (x2 - y2)^ (-1)) -4 * ((-1) + nu)* ((-3)* x3 - y3 + 2 *
        nu* (x3 + y3))* atan2(r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) +
        xlogy(x2 - y2, r1(y1, y2, y3) + x1 - y1) + xlogy((7 + 8 * ((-2) + nu)* nu)* (x2 +
          (-1)* y2), r2(y1, y2, y3) + x1 - y1))
    end

    function J2323(y1::R, y2::R, y3::R) where R
        @. (-1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (r1(y1, y2, y3) + r2(y1, y2, y3)^ (-1)* ((7 +
        8 * ((-2) + nu)* nu)* r2(y1, y2, y3)^ 2 + 2 * x3* y3) +2 * ((-3) + 4 * nu)* x3* acoth(
          r2(y1, y2, y3)^ (-1)* (x3 + y3)) + xlogy(2 * ((-3)* x3 - 2 * y3 + 6 * nu* (x3 + y3)
            -4 * nu^ 2 * (x3 + y3)), r2(y1, y2, y3) + x3 + y3))
    end

    function J3312(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x1 - y1)* (x2 - y2)* y3* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (
        -1)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (-1)* ((x1 - y1)^ 2 + (x2 + (
        -1)* y2)^ 2 + 2 * (x3 + y3)^ 2) -3 * x3* atan2(3 * x3, x1 - y1)
        -5 * x3* atan2(5 * x3, x2 - y2) + 12 * nu* x3* atan2((-3)* nu* x3, x2 + (
          -1)* y2) +4 * nu* x3* atan2(-nu* x3, x1 - y1) - 8 * nu^ 2 *
        x3* atan2(nu^ 2 * x3, x2 - y2) + 3 * y3* atan2((-3)* y3, x1 -
          y1) - 5 * y3* atan2(5 * y3, x2 - y2) + 12 * nu* y3* atan2((-3)*
            nu* y3, x2 - y2) - 4 * nu* y3* atan2(nu* y3, x1 - y1) - 8 *
        nu^ 2 * y3* atan2(nu^ 2 * y3, x2 - y2) + 2 * ((-1) + 2 * nu)* (x3 + (-1)
        * y3)* atan2(r1(y1, y2, y3)* (-x3 + y3), (x1 - y1)* (x2 - y2)) + 2 * (
        1 - 2 * nu)^ 2 * (x3 + y3)* atan2(r2(y1, y2, y3)* (x3 + y3), (x1 - y1)* (x2 + (-1)
        * y2)) +xlogy(-((-3) + 4 * nu)* (x1 - y1), r1(y1, y2, y3) + x2 - y2) +
        xlogy((5 + 4 * nu* ((-3) + 2 * nu))* (x1 - y1), r2(y1, y2, y3) + x2 - y2) +
        xlogy(-((-3) + 4 * nu)* (x2 - y2), r1(y1, y2, y3) + x1 - y1) + xlogy((5 +
          4 * nu* ((-3) + 2 * nu))* (x2 - y2), r2(y1, y2, y3) + x1 - y1))
    end

    function J3313(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x1 - y1)* y3* (x3 + y3)* ((x2 - y2)^ 2 + (x3 + y3)^ 2)^ (-1) + 5 *
        x2* atan2((-5)* x2, x1 - y1) - 3 * x2* atan2(3 * x2, x1 - y1)
        +4 * nu* x2* atan2(-nu* x2, x1 - y1) - 12 * nu* x2* atan2(
          3 * nu* x2, x1 - y1) + 8 * nu^ 2 * x2* atan2(-nu^ 2 * x2, x1 + (-1)
        * y1) - 4 * ((-1) + nu)* (x2 - y2)* atan2(r1(y1, y2, y3)* (x2 - y2), (x1 +
            (-1)* y1)* (x3 - y3)) -8 * ((-1) + nu)^ 2 * (x2 - y2)*
        atan2(r2(y1, y2, y3)* (x2 - y2), (x1 - y1)* (x3 + y3)) + 3 * y2* atan2((-3)
        * y2, x1 - y1) - 5 * y2* atan2(5 * y2, x1 - y1) + 12 * nu* y2*
        atan2((-3)* nu* y2, x1 - y1) - 4 * nu* y2* atan2(nu* y2, x1 + (-1)
        * y1) - 8 * nu^ 2 * y2* atan2(nu^ 2 * y2, x1 - y1) + xlogy((-4)*
          x3, r2(y1, y2, y3) - x1 + y1) + xlogy((-4)* ((-1) + nu)* (x1 - y1), r1(y1, y2, y3) + x3 + (-1)
        * y3) + xlogy((-8)* ((-1) + nu)^ 2 * (x1 - y1), r2(y1, y2, y3) + x3 + y3) + xlogy((-1)
        * ((-3) + 4 * nu)* (x3 - y3), r1(y1, y2, y3) + x1 - y1) + xlogy((-7)* x3
            -5 * y3 + 12 * nu* (x3 + y3) - 8 * nu^ 2 * (x3 + y3), r2(y1, y2, y3) + x1 - y1))
    end

    function J3323(y1::R, y2::R, y3::R) where R
        @. (1 / 16)* (1 - nu)^ (-1)* pi^ (-1)* G^ (-1)* (2 * r2(y1, y2, y3)^ (-1)* x3* (
        x2 - y2)* y3* (x3 + y3)* ((x1 - y1)^ 2 + (x3 + y3)^ 2)^ (-1) + 5 *
        x1* atan2((-5)* x1, x2 - y2) - 3 * x1* atan2(3 * x1, x2 - y2)
        +4 * nu* x1* atan2(-nu* x1, x2 - y2) - 12 * nu* x1* atan2(
          3 * nu* x1, x2 - y2) + 8 * nu^ 2 * x1* atan2(-nu^ 2 * x1, x2 + (-1)
        * y2) - 4 * ((-1) + nu)* (x1 - y1)* atan2(r1(y1, y2, y3)* (x1 - y1), (x2 +
            (-1)* y2)* (x3 - y3)) -8 * ((-1) + nu)^ 2 * (x1 - y1)*
        atan2(r2(y1, y2, y3)* (x1 - y1), (x2 - y2)* (x3 + y3)) + 3 * y1* atan2((-3)
        * y1, x2 - y2) - 5 * y1* atan2(5 * y1, x2 - y2) + 12 * nu* y1*
        atan2((-3)* nu* y1, x2 - y2) - 4 * nu* y1* atan2(nu* y1, x2 + (-1)
        * y2) - 8 * nu^ 2 * y1* atan2(nu^ 2 * y1, x2 - y2) + xlogy((-4)*
          x3, r2(y1, y2, y3) - x2 + y2) + xlogy((-4)* ((-1) + nu)* (x2 - y2), r1(y1, y2, y3) + x3 + (-1)
        * y3) + xlogy((-8)* ((-1) + nu)^ 2 * (x2 - y2), r2(y1, y2, y3) + x3 + y3) + xlogy((-1)
        * ((-3) + 4 * nu)* (x3 - y3), r1(y1, y2, y3) + x2 - y2) + xlogy((-7)* x3
            -5 * y3 + 12 * nu* (x3 + y3) - 8 * nu^ 2 * (x3 + y3), r2(y1, y2, y3) + x2 - y2))
    end

    function IU1(y1::R, y2::R, y3::R) where R
        @. (lambda * epsvkk + 2 * G * epsv11p) * J1123(y1, y2, y3)
        +2 * G * epsv12p * (J1223(y1, y2, y3) + J1113(y1, y2, y3))
        +2 * G * epsv13p * (J1323(y1, y2, y3) + J1112(y1, y2, y3))
        +(lambda * epsvkk + 2 * G * epsv22p) * J1213(y1, y2, y3)
        +2 * G * epsv23p * (J1212(y1, y2, y3) + J1313(y1, y2, y3))
        +(lambda * epsvkk + 2 * G * epsv33p) * J1312(y1, y2, y3)
    end

    function IU2(y1::R, y2::R, y3::R) where R
        @. (lambda * epsvkk + 2 * G * epsv11p) * J2123(y1, y2, y3)
        +2 * G * epsv12p * (J2223(y1, y2, y3) + J2113(y1, y2, y3))
        +2 * G * epsv13p * (J2323(y1, y2, y3) + J2112(y1, y2, y3))
        +(lambda * epsvkk + 2 * G * epsv22p) * J2213(y1, y2, y3)
        +2 * G * epsv23p * (J2212(y1, y2, y3) + J2313(y1, y2, y3))
        +(lambda * epsvkk + 2 * G * epsv33p) * J2312(y1, y2, y3)
    end

    function IU3(y1::R, y2::R, y3::R) where R
        @. (lambda * epsvkk + 2 * G * epsv11p) * J3123(y1, y2, y3)
        +2 * G * epsv12p * (J3223(y1, y2, y3) + J3113(y1, y2, y3))
        +2 * G * epsv13p * (J3323(y1, y2, y3) + J3112(y1, y2, y3))
        +(lambda * epsvkk + 2 * G * epsv22p) * J3213(y1, y2, y3)
        +2 * G * epsv23p * (J3212(y1, y2, y3) + J3313(y1, y2, y3))
        +(lambda * epsvkk + 2 * G * epsv33p) * J3312(y1, y2, y3)
    end

    function r1(y1::R, y2::R, y3::R) where R
        @. hypot(x1 - y1, x2 - y2, x3 - y3)
    end

    function r2(y1::R, y2::R, y3::R) where R
        @. hypot(x1 - y1, x2 - y2, x3 + y3)
    end

    lambda = G * 2 * nu / (1 - 2 * nu)
    epsvkk = epsv11p + epsv22p + epsv33p

    t1 = @. (x1 - q1) * cosd(theta) + (x2 - q2) * sind(theta)
    x2 = @. -(x1 - q1) * sind(theta) + (x2 - q2) * cosd(theta)
    x1 = t1

    u1 =
        IU1(L, T / 2, q3 + W) - IU1(L, -T / 2, q3 + W) + IU1(L, -T / 2, q3) - IU1(L, T / 2, q3)
        -IU1(0, T / 2, q3 + W) + IU1(0, -T / 2, q3 + W) - IU1(0, -T / 2, q3) + IU1(0, T / 2, q3)
    u2 =
        IU2(L, T / 2, q3 + W) - IU2(L, -T / 2, q3 + W) + IU2(L, -T / 2, q3) - IU2(L, T / 2, q3)
        -IU2(0, T / 2, q3 + W) + IU2(0, -T / 2, q3 + W) - IU2(0, -T / 2, q3) + IU2(0, T / 2, q3)
    u3 =
        IU3(L, T / 2, q3 + W) - IU3(L, -T / 2, q3 + W) + IU3(L, -T / 2, q3) - IU3(L, T / 2, q3)
        -IU3(0, T / 2, q3 + W) + IU3(0, -T / 2, q3 + W) - IU3(0, -T / 2, q3) + IU3(0, T / 2, q3)

    return cat(u1, u2, u3; dims=1)
end
