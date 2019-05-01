abstract type AbstracGreensFunction end
abstract type OkadaGreensFunction end

struct OkadaSelfieGreensFunction{T, F, A} <: OkadaGreensFunction
    tensor::T
    op::F
    alloc::A
end
