export max_velocity

"""
    max_velocity(t::AbstractVector, u::AbstractArray, getu::Function)

Return max velocity across the fault at each time step. A number of convenient interfaces for common
output are implemented.

## Arguments
- `t::AbstractVector`: vector of time steps
- `u::AbstractArray`: array of solution
- `getu::Function`: method for retrieving velocity section at each time step `i`, its signiture is `getu(u, i)`
"""
function max_velocity(t::AbstractVector, u::AbstractVector, getu::Function)
    x = similar(t)
    @fastmath @inbounds @simd for i = 1: length(x)
        x[i] = maximum(getu(u, i))
    end
    return x
end

function max_velocity(t::AbstractVector, u::AbstractVector{T}) where T<:ArrayPartition
    getu = (u, i) -> u[i].x[1] # assume velocity is stored in the first entry
    max_velocity(t, u, getu)
end

function max_velocity(sol::ODESolution)
    max_velocity(sol.t, sol.u)
end
