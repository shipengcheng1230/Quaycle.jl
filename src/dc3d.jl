function dc3d_fortran(x::T, y::T, z::T, α::T, dep::T, dip::T, al1::T, al2::T, aw1::T, aw2::T,
    disl1::T, disl2::T, disl3::T) where {T <: AbstractFloat}
    
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
    flag = Array{Int64}(1)
    
    # call okada's code which is renamed as "__dc3d__" (see binding rename in external/okada.f90)
    # input args tuple must be syntactically written instead of a variable assigned
    # macros could be used to simplify this in the future
    ccall((:__dc3d__, "./src/external/dc3d.so"), Void,
        (
            Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
            Ref{Int64},
        ),
        α, x, y, z, dep, dip, al1, al2, aw1, aw2, disl1, disl2, disl3,
        ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz,
        flag,
    )
    
    # results valid iff flag[1] == 0
    return (
        flag[1],
        ux[1], uy[1], uz[1],
        uxx[1], uyx[1], uzx[1],
        uxy[1], uyy[1], uzy[1],
        uxz[1], uyz[1], uzz[1]
    )
end


function dc3d_wrapper(
    coord::A, α::T, depth_fault::T, dip::T, range_strike::A, range_dip::A, dislocation::A,
    ) where {T <: AbstractFloat, A <: AbstractArray{T}}

    res = dc3d_fortran(coord..., α, depth_fault, dip, range_strike..., range_dip..., dislocation...)

    u = [res[2], res[3], res[4]]::Vector{T}
    ∇u = [
        res[5] res[8] res[11];
        res[6] res[9] res[12];
        res[7] res[10] res[13];
    ]::Array{T, 2}
    
    return (res[1], u, ∇u)
end
