module Base

abstract type PhysicalField end

struct FaultCell{TE, TA} <: PhysicalField where{TE <: AbstractFloat, TA <: AbstractArray{TE}}
    coord::TA
    u::TA
    ∇u::TA
end

end

struct Deformation
    u
    ∇u
end