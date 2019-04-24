## fault types

export DIPPING, STRIKING

abstract type AbstractFault end
abstract type PlaneFault <: AbstractFault end

struct DIPPING <: PlaneFault end
struct STRIKING <: PlaneFault end
