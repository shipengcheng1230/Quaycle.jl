"""
    max_velocity(sol)

Return max velocity across the fault at each time step.

## Arguments
- `sol`:: solution return by `DifferentialEquations.solve` given [`EarthquakeCycleProblem`](@ref).
"""
function max_velocity(sol)
    x = similar(sol.t)
    ndim = ndims(sol.u[1])

    for i = 1: length(x)
        x[i] = maximum(selectdim(sol.u[i], ndim, 1))
    end
    return x
end
