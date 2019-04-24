export max_velocity

"""
    max_velocity(t::AbstractVector, u::AbstractArray, getu::Function)

Return max velocity across the fault at each time step. A number of convenient interfaces for common
output are implemented.

## Arguments
- `t::AbstractVector`: vector of time steps
- `u::AbstractArray`: array of solution
- `getu::Function`: method for retrieving velocity section at each time step
"""
function max_velocity(t::AbstractVector, u::AbstractArray, getu::Function)
    x = similar(t)
    @fastmath @inbounds @simd for i = 1: length(x)
        x[i] = maximum(getu(u, i))
    end
    return x
end

function max_velocity(t::AbstractVector, u::AbstractArray{T, N}) where {T<:Number, N}
    @assert N == 3 || N == 4
    getu = N == 3 ? (u, i) -> view(u,:,1,i) : (u, i) -> view(u,:,:,1,i)
    max_velocity(t, u, getu)
end

function max_velocity(t::AbstractVector, u::AbstractArray{<:AbstractArray{T, N}}) where {T<:Number, N}
    getu = (u, i) -> selectdim(u[i], N, 1)
    max_velocity(t, u, getu)
end

function max_velocity(sol::ODESolution)
    max_velocity(sol.t, sol.u)
end
