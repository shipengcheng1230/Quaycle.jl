# Reference journal article:
# Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
# artefact-free solution.
# Submitted to Geophysical Journal International
#
# Copyright (c) 2014 Mehdi Nikkhoo
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the
# following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# I appreciate any comments or bug reports.
# 
# Mehdi Nikkhoo
# created: 2013.1.24
# Last modified: 2014.7.30
#
# VolcanoTectonics Research Group
# Section 2.1, Physics of Earthquakes and Volcanoes
# Department 2, Physics of the Earth
# Helmholtz Centre Potsdam
# German Research Centre for Geosciences (GFZ)
#
# email:
# mehdi.nikkhoo@gfz-potsdam.de
# mehdi.nikkhoo@gmail.com
# -------------------------------------------------------
# Translated by Pengcheng Shi (shipengcheng1230@gmail.com)

export td_disp_fs, td_disp_hs

const _ey = [0, 1, 0]
const _ez = [0, 0, 1]

function td_disp_hs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T) where {T, V}
    @assert (Z ≤ zero(T) && P1[3] ≤ zero(T) && P2[3] ≤ zero(T) && P3[3] ≤ zero(T))  "Half-space solution: Z coordinates must be negative!"
    A = transform_matrix(P1, P2, P3)
    ueMS, unMS, uvMS = _td_disp_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A')
    ueFSC, unFSC, uvFSC = _td_disp_harmonic_func(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A)
    P1[3] *= -one(T)
    P2[3] *= -one(T)
    P3[3] *= -one(T)
    allatsurface = P1[3] == zero(T) && P2[3] == zero(T) && P3[3] == zero(T)
    if !allatsurface
        A *= -one(T)
        A[3,1] *= -one(T)
        A[3,2] *= -one(T)
        A[3,3] *= -one(T)
    end
    ueIS, unIS, uvIS = _td_disp_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A')
    uvIS = ifelse(allatsurface, -uvIS, uvIS)
    ue = ueMS + ueIS + ueFSC
    un = unMS + unIS + unFSC
    uv = uvMS + uvIS + uvFSC
    allatsurface && begin ue *= -one(T); un *= -one(T); uv *= -one(T) end
    P1[3] *= -one(T)
    P2[3] *= -one(T)
    P3[3] *= -one(T)
    return ue, un, uv
end

function td_disp_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T) where {T, V}
    A = transform_matrix(P1, P2, P3)
    _td_disp_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A')
end

@inline function transform_matrix(P1::V, P2::V, P3::V) where V<:AbstractVector{T} where T
    Vnorm = cross(P2 - P1, P3 - P1)
    normalize!(Vnorm)
    Vstrike = cross(_ez, Vnorm)

    if norm(Vstrike) == zero(T)
        Vstrike = _ey * Vnorm[3]
        if P1[3] > zero(T)
            Vstrike *= -one(T)
        end
    end
    normalize!(Vstrike)
    Vdip = cross(Vnorm, Vstrike)
    hcat(Vnorm, Vstrike, Vdip)
end

@inline function _td_disp_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T, At::M) where {T, V, M}
    bx, by, bz = Ts, Ss, Ds
    p1, p2, p3 = [zeros(3) for _ in 1: 3]

    x, y, z = coord_trans(X - P2[1], Y - P2[2], Z - P2[3], At)
    p1[1], p1[2], p1[3] = coord_trans(P1[1] - P2[1], P1[2] - P2[2], P1[3] - P2[3], At)
    p3[1], p3[2], p3[3] = coord_trans(P3[1] - P2[1], P3[2] - P2[2], P3[3] - P2[3], At)

    e12 = normalize(p2 - p1)
    e13 = normalize(p3 - p1)
    e23 = normalize(p3 - p2)

    A = acos(dot(e12, e13))
    B = acos(dot(-e12, e23))
    C = acos(dot(e23, e13))

    Trimode = trimodefinder(y, z, x, p1, p2, p3)

    if Trimode == 1
        u1Tp, v1Tp, w1Tp = TDSetupD(x, y, z, A, bx, by, bz, nu, p1, -e13)
        u2Tp, v2Tp, w2Tp = TDSetupD(x, y, z, B, bx, by, bz, nu, p2, e12)
        u3Tp, v3Tp, w3Tp = TDSetupD(x, y, z, C, bx, by, bz, nu, p3, e23)
        u = u1Tp + u2Tp + u3Tp
        v = v1Tp + v2Tp + v3Tp
        w = w1Tp + w2Tp + w3Tp
    end

    if Trimode == -1
        u1Tn, v1Tn, w1Tn = TDSetupD(x, y, z, A, bx, by, bz, nu, p1, e13)
        u2Tn, v2Tn, w2Tn = TDSetupD(x, y, z, B, bx, by, bz, nu, p2, -e12)
        u3Tn, v3Tn, w3Tn = TDSetupD(x, y, z, C, bx, by, bz, nu, p3, -e23)
        u = u1Tn + u2Tn + u3Tn
        v = v1Tn + v2Tn + v3Tn
        w = w1Tn + w2Tn + w3Tn
    end

    if Trimode == 0
        u = NaN
        v = NaN
        w = NaN
    end

    a1, a2, a3 = -x, p1[2] - y, p1[3] - z
    b1, b2, b3 = -x, -y, -z
    c1, c2, c3 = -x, p3[2] - y, p3[3] - z
    na = hypot(a1, a2, a3)
    nb = hypot(b1, b2, b3)
    nc = hypot(c1, c2, c3)

    Fiy = a1 * (b2 * c3 - b3 * c2)- a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)
    Fix = na * nb * nc + (a1 * b1 + a2 * b2 + a3 * b3) * nc + (a1 * c1 + a2 * c2 + a3 * c3) * nb + (b1 * c1 + b2 * c2 + b3 * c3) * na
    Fi = -2 * atan(Fiy, Fix) / 4 / π
    u += bx * Fi
    v += by * Fi
    w += bz * Fi
    ue, un, uv = coord_trans(u, v, w, At')
end

@inline function _td_disp_harmonic_func(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T, A::M) where {T, V, M}
    bx, by, bz = Ts, Ss, Ds
    bX, bY, bZ = coord_trans(bx, by, bz, A)
    u1, v1, w1 = AngSetupFSC(X, Y, Z, bX, bY, bZ, P1, P2, nu)
    u2, v2, w2 = AngSetupFSC(X, Y, Z, bX, bY, bZ, P2, P3, nu)
    u3, v3, w3 = AngSetupFSC(X, Y, Z, bX, bY, bZ, P3, P1, nu)
    return u1 + u2 + u3, v1 + v2 + v3, w1 + w2 + w3
end

@inline function AngSetupFSC(X::T, Y::T, Z::T, bX::T, bY::T, bZ::T, PA::V, PB::V, nu::T) where {T, V}
    SideVec = PB - PA
    sv1, sv2, sv3 = PB[1] - PA[1], PB[2] - PA[2], PB[3] - PA[3]
    svr = hypot(sv1, sv2, sv3)

    beta = acos(-sv3 / svr)
    ey1 = [sv1, sv2, zero(T)]
    normalize!(ey1)
    ey3 = -_ez
    ey2 = cross(ey3, ey1)
    A = hcat(ey1, ey2, ey3)
    y1A, y2A, y3A = coord_trans(X - PA[1], Y - PA[2], Z - PA[3], A)
    y1AB, y2AB, y3AB = coord_trans(sv1, sv2, sv3, A)
    y1B = y1A - y1AB
    y2B = y2A - y2AB
    y3B = y3A - y3AB
    b1, b2, b3 = coord_trans(bX, bY, bZ, A)
    angle = ifelse(beta * y1A ≥ zero(T), -π + beta, beta)
    v1A, v2A, v3A = AngDisDispFSC(y1A, y2A, y3A, angle, b1, b2, b3, nu, -PA[3])
    v1B, v2B, v3B = AngDisDispFSC(y1B, y2B, y3B, angle, b1, b2, b3, nu, -PB[3])
    v1 = v1B - v1A
    v2 = v2B - v2A
    v3 = v3B - v3A
    ue, un, uv = coord_trans(v1, v2, v3, A')

    return ue, un, uv
end

@inline function AngDisDispFSC(y1::T, y2::T, y3::T, beta::T, b1::T, b2::T, b3::T, nu::T, a::T) where T
    sinB, cosB = sincos(beta)
    cotB = cot(beta)
    y3b = y3 + 2a
    z1b = y1 * cosB + y3b * sinB
    z3b = -y1 * sinB + y3b * cosB
    rb = hypot(y1, y2, y3b)
    Fib = 2 * atan(-y2 / (-(rb + y3b) * cot(beta/2) + y1))

    v1cb1 = (b1 / 4 / pi / (1 - nu) * (-2 * (1 - nu) * (1 - 2 * nu) * Fib * cotB ^ 2 + (1 - 2 * nu) * y2 /
        (rb + y3b) * ((1 - 2 * nu - a / rb) * cotB - y1 / (rb + y3b) * (nu + a / rb)) + (1 - 2 * nu) *
        y2 * cosB * cotB / (rb + z3b) * (cosB + a / rb) + a * y2 * (y3b - a) * cotB / rb ^ 3 + y2 *
        (y3b - a) / (rb * (rb + y3b)) * (-(1 - 2 * nu) * cotB + y1 / (rb + y3b) * (2 * nu + a / rb) +
        a * y1 / rb ^ 2) +y2 * (y3b - a) / (rb * (rb + z3b)) * (cosB / (rb + z3b) * ((rb *
        cosB + y3b) * ((1 - 2 * nu) * cosB - a / rb) * cotB + 2 * (1 - nu) * (rb * sinB - y1) * cosB) -
        a * y3b * cosB * cotB / rb ^ 2)))

    v2cb1 = ((b1 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * ((2 * (1 - nu) * cotB ^ 2 - nu) * log(rb + y3b) - (2 *
        (1 - nu) * cotB ^ 2 + 1 - 2 * nu)* cosB * log(rb + z3b)) -(1 - 2 * nu) / (rb + y3b) * (y1 *
        cotB * (1 - 2 * nu - a / rb) + nu * y3b - a + y2 ^ 2 / (rb + y3b) * (nu + a / rb)) -(1 - 2 *
        nu) * z1b * cotB / (rb + z3b) * (cosB + a / rb) - a * y1 * (y3b - a) * cotB / rb ^ 3 +
        (y3b - a) / (rb + y3b) * (-2 * nu + 1 / rb * ((1 - 2 * nu) * y1 * cotB - a) + y2 ^ 2 / (rb *
        (rb + y3b)) * (2 * nu + a / rb) + a * y2 ^ 2 / rb ^ 3) +(y3b - a) / (rb + z3b) * (cosB ^ 2 -
        1 / rb * ((1 - 2 * nu) * z1b * cotB + a * cosB) + a * y3b * z1b * cotB / rb ^ 3 - 1 / (rb *
        (rb + z3b)) * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB / rb * (rb * cosB + y3b))))))

    v3cb1 = (b1 / 4 / pi / (1 - nu) * (2 * (1 - nu) * (((1 - 2 * nu) * Fib * cotB) + (y2 / (rb + y3b) * (2 *
        nu + a / rb)) -(y2 * cosB / (rb + z3b) * (cosB + a / rb))) +y2 * (y3b - a) / rb * (2 *
        nu / (rb + y3b) + a / rb ^ 2) +y2 * (y3b - a) * cosB / (rb * (rb + z3b)) * (1 - 2 * nu -
        (rb * cosB + y3b) / (rb + z3b) * (cosB + a / rb) - a * y3b / rb ^ 2)))

    v1cb2 = (b2 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * ((2 * (1 - nu) * cotB ^ 2 + nu) * log(rb + y3b) - (2 *
        (1 - nu) * cotB ^ 2 + 1)* cosB * log(rb + z3b)) +(1 - 2 * nu) / (rb + y3b) * (-(1 - 2 * nu) *
        y1 * cotB + nu * y3b - a + a * y1 * cotB / rb + y1 ^ 2 / (rb + y3b) * (nu + a / rb)) -(1 - 2 *
        nu)* cotB / (rb + z3b) * (z1b * cosB - a * (rb * sinB - y1) / (rb * cosB)) - a * y1 *
        (y3b - a) * cotB / rb ^ 3 + (y3b - a) / (rb + y3b) * (2 * nu + 1 / rb * ((1 - 2 * nu) * y1 *
        cotB + a) -y1 ^ 2 / (rb * (rb + y3b)) * (2 * nu + a / rb) - a * y1 ^ 2 / rb ^ 3) +(y3b - a) *
        cotB / (rb + z3b) * (-cosB * sinB + a * y1 * y3b / (rb ^ 3 * cosB) + (rb * sinB - y1) /
        rb * (2 * (1 - nu) * cosB - (rb * cosB + y3b) / (rb + z3b) * (1 + a / (rb * cosB))))))

    v2cb2 = (b2 / 4 / pi / (1 - nu) * (2 * (1 - nu) * (1 - 2 * nu) * Fib * cotB ^ 2 + (1 - 2 * nu) * y2 /
        (rb + y3b) * (-(1 - 2 * nu - a / rb) * cotB + y1 / (rb + y3b) * (nu + a / rb)) - (1 - 2 * nu) *
        y2 * cotB / (rb + z3b) * (1 + a / (rb * cosB)) - a * y2 * (y3b - a) * cotB / rb ^ 3 + y2 *
        (y3b - a) / (rb * (rb + y3b)) * ((1 - 2 * nu) * cotB - 2 * nu * y1 / (rb + y3b) - a * y1 / rb *
        (1 / rb + 1 / (rb + y3b))) +y2 * (y3b - a) * cotB / (rb * (rb + z3b)) * (-2 * (1 - nu) *
        cosB + (rb * cosB + y3b) / (rb + z3b) * (1 + a / (rb * cosB)) + a * y3b / (rb ^ 2 * cosB))))

    v3cb2 = (b2 / 4 / pi / (1 - nu) * (-2 * (1 - nu) * (1 - 2 * nu) * cotB * (log(rb + y3b) - cosB *
        log(rb + z3b)) -2 * (1 - nu) * y1 / (rb + y3b) * (2 * nu + a / rb) + 2 * (1 - nu) * z1b / (rb +
        z3b) * (cosB + a / rb) + (y3b - a) / rb * ((1 - 2 * nu) * cotB - 2 * nu * y1 / (rb + y3b) - a *
        y1 / rb ^ 2) -(y3b - a) / (rb + z3b) * (cosB * sinB + (rb * cosB + y3b) * cotB / rb *
        (2 * (1 - nu) * cosB - (rb * cosB + y3b) / (rb + z3b)) + a / rb * (sinB - y3b * z1b /
        rb ^ 2 - z1b * (rb * cosB + y3b) / (rb * (rb + z3b))))))

    v1cb3 = (b3 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * (y2 / (rb + y3b) * (1 + a / rb) - y2 * cosB / (rb +
        z3b) * (cosB + a / rb)) -y2 * (y3b - a) / rb * (a / rb ^ 2 + 1 / (rb + y3b)) + y2 *
        (y3b - a) * cosB / (rb * (rb + z3b)) * ((rb * cosB + y3b) / (rb + z3b) * (cosB + a /
        rb) +a * y3b / rb ^ 2)))

    v2cb3 = (b3 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * (-sinB * log(rb + z3b) - y1 / (rb + y3b) * (1 + a /
        rb) +z1b / (rb + z3b) * (cosB + a / rb)) +y1 * (y3b - a) / rb * (a / rb ^ 2 + 1 / (rb +
        y3b)) -(y3b - a) / (rb + z3b) * (sinB * (cosB - a / rb) + z1b / rb * (1 + a * y3b /
        rb ^ 2) -1 / (rb * (rb + z3b)) * (y2 ^ 2 * cosB * sinB - a * z1b / rb * (rb * cosB + y3b)))))

    v3cb3 = (b3 / 4 / pi / (1 - nu) * (2 * (1 - nu) * Fib + 2 * (1 - nu) * (y2 * sinB / (rb + z3b) * (cosB +
        a / rb)) +y2 * (y3b - a) * sinB / (rb * (rb + z3b)) * (1 + (rb * cosB + y3b) / (rb +
        z3b) * (cosB + a / rb) + a * y3b / rb ^ 2)))

    return v1cb1 + v1cb2 + v1cb3, v2cb1 + v2cb2 + v2cb3, v3cb1 + v3cb2 + v3cb3
end

@inline function coord_trans(x1::T, x2::T, x3::T, A::M) where {T<:Number, M<:AbstractMatrix}
    r1 = A[1,1] * x1 + A[1,2] * x2 + A[1,3] * x3
    r2 = A[2,1] * x1 + A[2,2] * x2 + A[2,3] * x3
    r3 = A[3,1] * x1 + A[3,2] * x2 + A[3,3] * x3
    r1, r2, r3
end

@inline function trimodefinder(x::T, y::T, z::T, p1::V, p2::V, p3::V) where {T, V}
    a = ((p2[3] - p3[3]) * (x - p3[2]) + (p3[2] - p2[2]) * (y - p3[3])) /
        ((p2[3] - p3[3]) * (p1[2] - p3[2]) + (p3[2] - p2[2]) * (p1[3] - p3[3]))
    b = ((p3[3] - p1[3]) * (x - p3[2]) + (p1[2] - p3[2]) * (y - p3[3])) /
        ((p2[3] - p3[3]) * (p1[2] - p3[2]) + (p3[2] - p2[2]) * (p1[3] - p3[3]))
    c = 1 - a - b
    trimode = 1
    if ((a ≤ 0) && (b > c) && (c > a)) || ((b ≤ 0) && (c > a) && (a > b)) || ((c ≤ 0) && (a > b) && (b > c))
        trimode = -1
    end
    if ((a == 0) && (b ≥ 0) && (c ≥ 0)) || ((a ≥ 0) && (b == 0) && (c ≥ 0)) || ((a ≥ 0) && (b ≥ 0) && (c == 0))
        trimode = 0
    end
    if (trimode == 0) && (z ≠ 0)
        trimode = 1
    end
    return trimode
end

@inline function TDSetupD(x::T, y::T, z::T, alpha::T, bx::T, by::T, bz::T, nu::T, TriVertex::V, SideVec::V) where {T, V}
    y1 = SideVec[3] * (y - TriVertex[2]) - SideVec[2] * (z - TriVertex[3])
    z1 = SideVec[2] * (y - TriVertex[2]) + SideVec[3] * (z - TriVertex[3])
    by1 = SideVec[3] * by - SideVec[2] * bz
    bz1 = SideVec[2] * by + SideVec[3] * bz
    u, v0, w0 = AngDisDisp(x, y1, z1, -π + alpha, bx, by1, bz1, nu)
    v = SideVec[3] * v0 + SideVec[2] * w0
    w = -SideVec[2] * v0 + SideVec[3] * w0
    return u, v, w
end

@inline function AngDisDisp(x::T, y::T, z::T, alpha::T, bx::T, by::T, bz::T, nu::T) where T
    sinA, cosA = sincos(alpha)
    eta = y * cosA - z * sinA
    zeta = y * sinA + z * cosA
    r = hypot(x, y, z)

    zeta = ifelse(zeta > r, r, zeta)
    z = ifelse(z > r, r, z)

    ux = bx / 8 / pi / (1 - nu) * (x * y / r / (r - z) - x * eta / r / (r - zeta))
    vx = bx / 8 / pi / (1 - nu) * (eta * sinA / (r - zeta) - y * eta / r / (r - zeta) + y ^ 2 / r / (r - z) + (1 - 2 * nu) * (cosA * log(r - zeta) - log(r - z)))
    wx = bx / 8 / pi / (1 - nu) * (eta * cosA / (r - zeta) - y / r - eta * z / r / (r - zeta) - (1 - 2 * nu) * sinA * log(r - zeta))

    uy = by / 8 / pi / (1 - nu) * (x ^ 2 * cosA / r / (r - zeta) - x ^ 2 / r / (r - z) - (1 - 2 * nu) * (cosA * log(r - zeta) - log(r - z)))
    vy = by * x / 8 / pi / (1 - nu) * (y * cosA / r / (r - zeta) - sinA * cosA / (r - zeta) - y / r / (r - z))
    wy = by * x / 8 / pi / (1 - nu) * (z * cosA / r / (r - zeta) - cosA ^ 2 / (r - zeta) + 1 / r)

    uz = bz * sinA / 8 / pi / (1 - nu) * ((1 - 2 * nu) * log(r - zeta) - x ^ 2 / r / (r - zeta))
    vz = bz * x * sinA / 8 / pi / (1 - nu) * (sinA / (r - zeta) - y / r / (r - zeta))
    wz = bz * x * sinA / 8 / pi / (1 - nu) * (cosA / (r - zeta) - z / r / (r - zeta))

    return ux + uy + uz, vx + vy + vz, wx + wy + wz
end
