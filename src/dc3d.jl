function dc3d(x::T, y::T, z::T, α::T, dep::T, dip::T, al1::T, al2::T, aw1::T, aw2::T,
    disl1::T, disl2::T, disl3::T) where {T<:AbstractFloat}

    # initial return values
    ux = Array{Float64}(1)
    uy = Array{Float64}(1)
    uz = Array{Float64}(1)
    uxx = Array{Float64}(1)
    uyx = Array{Float64}(1)
    uzx = Array{Float64}(1)
    uxy = Array{Float64}(1)
    uyy = Array{Float64}(1)
    uzy = Array{Float64}(1)
    uxz = Array{Float64}(1)
    uyz = Array{Float64}(1)
    uzz = Array{Float64}(1)
    iret = Array{Int64}(1)

    # call okada's code
    # input args tuple must be syntactically written instead of a variable assigned
    # macros could be used to simplify this
    ccall((:__dc3d__, "./src/external/dc3d.so"), Void,
        (Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
         Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
         Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
         Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
         Ref{Int64}),
        α, x, y, z, dep, dip, al1, al2, aw1, aw2, disl1, disl2, disl3,
        ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz,
        iret
    )

    # results valid iif iret[1] == 0
    if iret[1] == 0
        return (ux[1], uy[1], uz[1],
                uxx[1], uyx[1], uzx[1],
                uxy[1], uyy[1], uzy[1],
                uxz[1], uyz[1], uzz[1])
    else
        # use Value Type for type stability
        return ntuple(_ -> NaN, Val{12})
    end
end
