module JuEQ

include("rsf.jl")
include("fault.jl")
include("bem.jl")
include("utils.jl")

export NormalFault, ThrustFault, StrikeSlipFault
export PlaneFaultDomain
export fault

export DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export CForm, RForm
export friction
export MaterialProperties

export discretize, properties
export EarthquakeCycleProblem
export dc3d_okada

export max_velocity

end
