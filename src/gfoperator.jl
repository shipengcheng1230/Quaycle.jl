## operators of greens function

abstract type ODEStateVariable{dim} end

struct ODEState_1D{T<:AbstractVecOrMat{<:Real}} <: ODEStateVariable{1}
    dτ_dt::T
    relv::T
end

struct ODEState_2D{T<:AbstractArray{<:Real}, U<:AbstractArray{<:Complex}, P<:Plan} <: ODEStateVariable{2}
    dτ_dt::T
    relv::T
    dτ_dt_dft::U
    relv_dft::U
    dτ_dt_buffer::T
    pf::P
end

create_ode_state_vars(nξ::Integer; T1=Float64) = ODEState_1D([Vector{T1}(undef, nξ) for _ in 1: 2]...)

function create_ode_state_vars(nx::I, nξ::I; T1=Float64) where {I <: Integer}
    FFTW.set_num_threads(FFT_CONFIGS["FFT_NUM_THREADS"])
    x1 = Matrix{T1}(undef, 2 * nx - 1, nξ)
    p1 = plan_rfft(x1, 1, flags=FFT_CONFIGS["FFT_FLAG"])

    tvar = ODEState_2D(
        Matrix{T1}(undef, nx, nξ),
        zeros(T1, 2nx-1, nξ), # for relative velocity
        [Matrix{Complex{T1}}(undef, nx, nξ) for _ in 1: 2]...,
        Matrix{T1}(undef, 2nx-1, nξ),
        p1)
end

create_ode_state_vars(gd::BoundaryElementGrid{1}) = create_ode_state_vars(gd.nξ; T1=typeof(gd.Δξ))
create_ode_state_vars(gd::BoundaryElementGrid{2}) = create_ode_state_vars(gd.nx, gd.nξ; T1=typeof(gd.Δx))

function derivations!(du, u, mp::PlaneMaterialProperties{dim}, tvar::ODEStateVariable{dim}, se::StateEvolutionLaw, fform::FrictionLawForm) where {dim}
    _ndim = dim + 1
    v = selectdim(u, _ndim, 1)
    θ = selectdim(u, _ndim, 2)
    dv = selectdim(du, _ndim, 1)
    dθ = selectdim(du, _ndim, 2)

    variable_trim!(v, θ)
    dτ_dt!(mp, tvar, v)
    dv_dθ_dt!(fform, dv, dθ, v, θ, mp, tvar, se)
end

@inline function variable_trim!(v::T, θ::T) where {T<:AbstractVecOrMat{<:Float64}}
    clamp!(θ, 0., Inf)
end

function gf_operator_kernel(::Val{:okada}, gf::AbstractArray{T, 2}, v::AbstractVector, v0::T) where T

end


function gf_operator_kernel(::Val{:okada}, gf::AbstractArray{T, 3}, v::AbstractMatrix, v0::T) where T

end


@inline function dτ_dt!(mp::PlaneMaterialProperties{1}, tvar::ODEStateVariable{1}, v::AbstractArray)
    @fastmath @inbounds @simd for i = 1: prod(mp.dims)
        tvar.relv[i] = mp.vpl - v[i]
    end
    mul!(tvar.dτ_dt, mp.k, tvar.relv)
end

@inline function dτ_dt!(mp::PlaneMaterialProperties{2}, tvar::ODEStateVariable{2}, v::AbstractArray{T}) where {T<:Number}
    @fastmath @inbounds for j = 1: mp.dims[2]
        @simd for i = 1: mp.dims[1]
            tvar.relv[i,j] = mp.vpl - v[i,j]
        end
    end
    mul!(tvar.relv_dft, tvar.pf, tvar.relv)
    fill!(tvar.dτ_dt_dft, zero(T))

    @fastmath @inbounds for l = 1: mp.dims[2]
        for j = 1: mp.dims[2]
            @simd for i = 1: mp.dims[1]
                tvar.dτ_dt_dft[i,j] += mp.k[i,j,l] * tvar.relv_dft[i,l]
            end
        end
    end

    ldiv!(tvar.dτ_dt_buffer, tvar.pf, tvar.dτ_dt_dft)

    @fastmath @inbounds for j = 1: mp.dims[2]
        @simd for i = 1: mp.dims[1]
            tvar.dτ_dt[i,j] = tvar.dτ_dt_buffer[i,j]
        end
    end
end
