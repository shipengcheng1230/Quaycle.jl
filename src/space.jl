export compose
export DIPPING, STRIKING

abstract type AbstractFaultType end
abstract type PlaneFault <: AbstractFaultType end
abstract type FlatPlaneFault <: PlaneFault end

"Dipping, indicate dislocation occurs at downdip direction."
struct DIPPING <: FlatPlaneFault end
"Striking, indicate dislocation occurs at strike direction."
struct STRIKING <: FlatPlaneFault end

"Compose multiple slip vectors on the fault plane."
struct CompositePlaneFault{V} <: PlaneFault
    components::V
end
