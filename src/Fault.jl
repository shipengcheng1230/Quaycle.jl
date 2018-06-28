module Fault

abstract type FaultDomain end
abstract type Fault2D <: FaultDomain end
abstract type TransformFault <: Fault2D end

struct UniformTransformFault{T, A} <: TransformFault where {T <: Number, A <: AbstractVector{T}}
    xi::A
    zi::A
end

end
