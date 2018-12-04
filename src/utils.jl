"""
    max_velocity(t::AbstractVector, u::AbstractArray, getu::Function)

Return max velocity across the fault at each time step. A number of convenient interfaces for common
output are implemented.

## Arguments
- t::AbstractVector: vector of time steps
- u::AbstractArray: array of solution
- getu::Function: method for retrieving velocity section at each time step
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

"Calculate moment magnitude."
moment_magnitude(μ::T, d::T, A::T) where {T<:Number} = μ * d * A

"""
    DECallbackSaveToFile(iot::IOStream, iou::IOStream)

Construct a functional callback to write `ODESolution` (`t` & `u`) into file. The reason to separate `t` and `u` is
for more easily reshape `u` w.r.t grids specification.

**Note**
It is strongly not recommended to use "skipping" scheme (by defining `thrd` and `dts(a)` for each case) when solution is too
oscillated.
"""
function DECallbackSaveToFile(iot::IOStream, iou::IOStream)

    function __save(u, t, integrator)
        write(iot, t)
        write(iou, u)
    end
    FunctionCallingCallback(__save; func_everystep=true)
end

function DECallbackSaveToFile(iot::IOStream, iou::IOStream, thrd::Number; dts=0.1, dta=3.1536e6)

    function __save(u::AbstractArray{T, N}, t, integrator) where {T, N}
        maxv = maximum(selectdim(u, N, 1))
        maxvprev = maximum(selectdim(integrator.uprev, N, 1))
        if (maxv ≥ thrd && t - integrator.tprev ≥ dts) || (maxvprev < thrd && t - integrator.tprev ≥ dta)
            write(iot, t)
            write(iou, u)
        end
    end
    FunctionCallingCallback(__save; func_everystep=true)
end
