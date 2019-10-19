@inline function coord_trans(x1::T, x2::T, x3::T, A::M) where {T<:Number, M<:AbstractMatrix}
    r1 = A[1,1] * x1 + A[1,2] * x2 + A[1,3] * x3
    r2 = A[2,1] * x1 + A[2,2] * x2 + A[2,3] * x3
    r3 = A[3,1] * x1 + A[3,2] * x2 + A[3,3] * x3
    r1, r2, r3
end

@inline function coord_trans(x1::T, x2::T, x3::T, A::M) where {T<:AbstractVector, M<:AbstractMatrix}
    r1, r2, r3 = similar(x1), similar(x2), similar(x3)
    @inbounds @simd for i ∈ eachindex(r1)
        r1[i] = A[1,1] * x1[i] + A[1,2] * x2[i] + A[1,3] * x3[i]
        r2[i] = A[2,1] * x1[i] + A[2,2] * x2[i] + A[2,3] * x3[i]
        r3[i] = A[3,1] * x1[i] + A[3,2] * x2[i] + A[3,3] * x3[i]
    end
    r1, r2, r3
end


@inline function AngDisDisp(x::U, y::U, z::U, alpha::T, bx::T, by::T, bz::T, nu::T) where {T<:Number, U}
    sinA, cosA = sincos(alpha)
    eta = @. y * cosA - z * sinA
    zeta = @. y * sinA + z * cosA
    r = hypot.(x, y, z)

    if isa(U, Number)
        zeta = ifelse(zeta > r, r, zeta)
        z = ifelse(z > r, r, z)
    elseif isa(U, AbstractArray)
        bool = zeta .> r
        zeta[bool] .= r[bool]
        @. bool = z > r
        z[bool] .= r[bool]
    end

    ux = @. bx / 8 / pi / (1 - nu) * (x * y / r / (r - z) - x * eta / r / (r - zeta))
    vx = @. bx / 8 / pi / (1 - nu) * (eta * sinA / (r - zeta) - y * eta / r / (r - zeta) + y ^ 2 / r / (r - z) + (1 - 2 * nu) * (cosA * log(r - zeta) - log(r - z)))
    wx = @. bx / 8 / pi / (1 - nu) * (eta * cosA / (r - zeta) - y / r - eta * z / r / (r - zeta) - (1 - 2 * nu) * sinA * log(r - zeta))

    uy = @. by / 8 / pi / (1 - nu) * (x ^ 2 * cosA / r / (r - zeta) - x ^ 2 / r / (r - z) - (1 - 2 * nu) * (cosA * log(r - zeta) - log(r - z)))
    vy = @. by * x / 8 / pi / (1 - nu) * (y * cosA / r / (r - zeta) - sinA * cosA / (r - zeta) - y / r / (r - z))
    wy = @. by * x / 8 / pi / (1 - nu) * (z * cosA / r / (r - zeta) - cosA ^ 2 / (r - zeta) + 1 / r)

    uz = @. bz * sinA / 8 / pi / (1 - nu) * ((1 - 2 * nu) * log(r - zeta) - x ^ 2 / r / (r - zeta))
    vz = @. bz * x * sinA / 8 / pi / (1 - nu) * (sinA / (r - zeta) - y / r / (r - zeta))
    wz = @. bz * x * sinA / 8 / pi / (1 - nu) * (cosA / (r - zeta) - z / r / (r - zeta))

    return ux + uy + uz, vx + vy + vz, wx + wy + wz
end
