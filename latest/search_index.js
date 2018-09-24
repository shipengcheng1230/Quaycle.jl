var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#JuEQ.jl-Documentation-1",
    "page": "Home",
    "title": "JuEQ.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "This is a suite for numerically simulating earthquake sequences in Julia. The purpose of this package is to provide efficient Julia implementations for simulations in the field of earthquake physics.Features of this package are listed as below:Rate-State Friction Law\nOkada\'s Dislocation Method\nBoundary Element Method (Quasi-dynamic)Features to be implemented:Viscoelastic Relaxation (priority)\nFully Elastodynamic Effect\nOff-Fault Materials Effect"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Get the latest version with Julia\'s package manager:] add https://github.com/shipengcheng1230/JuEQ.jlTo load the package:using JuEQ"
},

{
    "location": "index.html#Acknowledgements-1",
    "page": "Home",
    "title": "Acknowledgements",
    "category": "section",
    "text": "The simulation of episodic seismic and slow slip events using boundary-element-method is largely benifit from Yajing Liu original Fortran code."
},

{
    "location": "quasi_dynamic_intro.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "quasi_dynamic_intro.html#Quasi-dynamic-Simulation-1",
    "page": "Introduction",
    "title": "Quasi-dynamic Simulation",
    "category": "section",
    "text": ""
},

{
    "location": "quasi_dynamic_intro.html#Basic-Theory-1",
    "page": "Introduction",
    "title": "Basic Theory",
    "category": "section",
    "text": "The governing equation is that at every time step, shear stress across the fault plane equals to frictional force plus a radiation damping term for approximating wave propagation effect:τ = σf + ηVHere μ is shear stress across the fault plain. Using Okada\'s dislocation theory, it can be shown as:τ = K  δwhere K is the so-called stiffness tensor, depicting relationship between displacements at one position regarding to dislocations somewhere else. δ is the dislocation, i.e. displacement at everywhere on the fault.  denotes tensor contraction.Back to f, we use rate-and-state frictional law to calculate its value, specifically as below:f(V θ) = f_0 + a lnfracVV_0 + b lnleft(fracV_0 θLright)where f_0 and V_0 are reference friction coefficient and velocity, V and θ are velocity and state variable based on which frictional force is. a and b are two frictional parameters denoting contributions each of which comes from velocity and state variable respectively.Sometimes people use regularized form to avoid infinity when V  0, namely:f(V θ) = a sinh ^-1leftfracV2V_0 expfracf_0 + b lnleft(V_0 θLright)arightThere are many state evolution law that describes how state variable θ changes with time, one of which that most widely used is Dieterich law:fracmathrmdθmathrmdt = 1 - fracV θLFurther, η is a damping coefficient whose value is often chosen as μ  2mathrmVs where μ is shear modulus and mathrmVs shear wave velocity. σ is the effective normal stress.To simulate how fault evolves with time, we then take the derivative of the governing equation:fracmathrmd τmathrmd t = fracmathrmd f(V θ)mathrmd t + η fracmathrmd Vmathrmd t\n= fracmathrmd fmathrmd V fracmathrmd Vmathrmd t + fracmathrmd fmathrmd θ fracmathrmd θmathrmd t + η fracmathrmd Vmathrmd tThus we arrive at:fracmathrmd Vmathrmd t = fracfracmathrmd τmathrmd t - fracmathrmd τmathrmd θ fracmathrmd θmathrmd tfracmathrmd τmathrmd V + ηwhere fracmathrmd τmathrmd t = K  (mathrmV_pl - V) where mathrmV_pl is the plate rate.Hence, with both derivatives of velocity V and state variable θ, we are able to discover how fault evolves with various parameters settings."
},

{
    "location": "quasi_dynamic_example.html#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "quasi_dynamic_example.html#Quasi-dynamic-Simulation-1",
    "page": "Examples",
    "title": "Quasi-dynamic Simulation",
    "category": "section",
    "text": ""
},

{
    "location": "quasi_dynamic_example.html#Example-1:-1D-fault-1",
    "page": "Examples",
    "title": "Example 1: 1D fault",
    "category": "section",
    "text": "note: Note\nThis example is from Benchmark Problem 1(hence referred as BP1)"
},

{
    "location": "quasi_dynamic_example.html#Parameters-Setting-1",
    "page": "Examples",
    "title": "Parameters Setting",
    "category": "section",
    "text": ""
},

{
    "location": "public_interface.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "public_interface.html#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": ""
},

{
    "location": "public_interface.html#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public_interface.md\"]"
},

{
    "location": "public_interface.html#JuEQ.DieterichStateLaw",
    "page": "Public",
    "title": "JuEQ.DieterichStateLaw",
    "category": "type",
    "text": "\\frac{\\mathrm{d}θ}{\\mathrm{d}t} = 1 - \\frac{v θ}{L}`\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.PrzStateLaw",
    "page": "Public",
    "title": "JuEQ.PrzStateLaw",
    "category": "type",
    "text": "fracmathrmdθmathrmdt = 1 - (fracv θ2L)^2\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.RuinaStateLaw",
    "page": "Public",
    "title": "JuEQ.RuinaStateLaw",
    "category": "type",
    "text": "fracmathrmdθmathrmdt = -fracv θL * logfracv θL\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.EarthquakeCycleProblem-Union{Tuple{dim}, Tuple{BoundaryElementGrid,PlaneMaterialProperties{dim,U,P,T} where T<:Number where P<:AbstractArray where U<:(Union{AbstractArray{T,1}, AbstractArray{T,2}} where T),Any,Any}} where dim",
    "page": "Public",
    "title": "JuEQ.EarthquakeCycleProblem",
    "category": "method",
    "text": "EarthquakeCycleProblem(gd::BoundaryElementGrid, p::PlaneMaterialProperties, u0, tspan; se=DieterichStateLaw(), fform=CForm())\n\nReturn an ODEProblem that encapsulate all the parameters and functions required for simulation. For the entailing usage, please refer DifferentialEquations.jl\n\nArguments\n\ngd::BoundaryElementGrid: grid set by user given a fault domain.\np::PlaneMaterialProperties: material profile w.r.t. the fault grid.\nu0: initial condition, should be organized such that the first of last dim is velocity while the 2nd of last dim is state.\ntspan: time interval to be simulated.\nse::StateEvolutionLaw: state evolution law to be applied.\nfform::FrictionLawForm: forms of frictional law to be applied.\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.dc3d_okada-Union{Tuple{A}, Tuple{T}, Tuple{T,T,T,T,T,T,Union{SubArray, A},Union{SubArray, A},A}} where A<:Union{AbstractArray{T,1}, AbstractArray{T,2}} where T<:Number",
    "page": "Public",
    "title": "JuEQ.dc3d_okada",
    "category": "method",
    "text": "Calculate displacements and gradient of displacements due to a dislocation in an elastic isotropic halfspace. See dc3d for details.\n\ntest/test_okada.dat is obtained using DC3dfortran\n\nAn example wrapper for DC3D in julia as below:\n\nfunction dc3d_fortran(x::T, y::T, z::T, α::T, dep::T, dip::T, al1::T, al2::T, aw1::T, aw2::T,\n    disl1::T, disl2::T, disl3::T) where {T <: AbstractFloat}\n\n    # initial return values\n    # `RefValue{T}` may be also viable other than `Array{T, 1}`\n    ux = Array{Float64}(1)\n    uy = Array{Float64}(1)\n    uz = Array{Float64}(1)\n    uxx = Array{Float64}(1)\n    uyx = Array{Float64}(1)\n    uzx = Array{Float64}(1)\n    uxy = Array{Float64}(1)\n    uyy = Array{Float64}(1)\n    uzy = Array{Float64}(1)\n    uxz = Array{Float64}(1)\n    uyz = Array{Float64}(1)\n    uzz = Array{Float64}(1)\n    iret = Array{Int64}(1)\n\n    # call okada\'s code which is renamed as \"__dc3d__\" (see binding rename shown below)\n    # input args tuple must be syntactically written instead of a variable assigned\n    # macros could be used to simplify this in the future\n    ccall((:__dc3d__, \"dc3d.so\"), Void,\n        (\n            Ref{Float64},\n            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},\n            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},\n            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ref{Int64},\n        ),\n        α, x, y, z, dep, dip, al1, al2, aw1, aw2, disl1, disl2, disl3,\n        ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz,\n        iret,\n    )\n\n    # results valid iff iret[1] == 0\n    return (\n        iret[1],\n        ux[1], uy[1], uz[1],\n        uxx[1], uyx[1], uzx[1],\n        uxy[1], uyy[1], uzy[1],\n        uxz[1], uyz[1], uzz[1]\n    )\nend\n\nThe corresponding fortran module is:\n\nMODULE okada\n  USE, INTRINSIC :: iso_c_binding\n  IMPLICIT NONE\nCONTAINS\n\n  SUBROUTINE dc3d_wrapper(&\n       & alpha, &\n       & x, y, z, &\n       & depth, dip, &\n       & al1, al2, &\n       & aw1, aw2, &\n       & disl1, disl2, disl3, &\n       & ux, uy, uz, &\n       & uxx, uyx, uzx, &\n       & uxy, uyy, uzy, &\n       & uxz, uyz, uzz, &\n       & iret) BIND(C, NAME=\'__dc3d__\')\n\n    REAL*8 :: &\n         & alpha, &\n         & x, y, z, &\n         & depth, dip, &\n         & al1, al2, &\n         & aw1, aw2, &\n         & disl1, disl2, disl3, &\n         & ux, uy, uz, &\n         & uxx, uyx, uzx, &\n         & uxy, uyy, uzy, &\n         & uxz, uyz, uzz\n\n    INTEGER*8 :: iret\n\n    CALL dc3d(&\n         & alpha, &\n         & x, y, z, &\n         & depth, dip, &\n         & al1, al2, &\n         & aw1, aw2, &\n         & disl1, disl2, disl3, &\n         & ux, uy, uz, &\n         & uxx, uyx, uzx, &\n         & uxy, uyy, uzy, &\n         & uxz, uyz, uzz, &\n         & iret)\n\n  END SUBROUTINE dc3d_wrapper\n\nEND MODULE okada\n\nA sample of makefile is as below:\n\n# Build Okada\'s code for calculating deformation due to a fault model\n#\nCC = gfortran\nCFLAGS = -fPIC -w -O3\nLDFLAGS = -shared\n\nSRCS = dc3d.f okada.f90\nOBJS = $(SRCS:.c=.o)\n\nTARGET = dc3d.so\n\n$(TARGET): $(OBJS)\n    $(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS)\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.discretize-Union{Tuple{T}, Tuple{ftype}, Tuple{PlaneFaultDomain{ftype,1,T} where T,T}} where T where ftype",
    "page": "Public",
    "title": "JuEQ.discretize",
    "category": "method",
    "text": "discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ::T; ax_ratio=50*one(T)::T)\n\nGenerate the grid for given 1D fault domain.\n\nArguments\n\nΔξ: grid space along-downdip\nax_ratio::Number: ration of along-strike length agsinst along-downdip length for mimicing an extended   2d (x & ξ) fault represented by 1d (ξ) domain. Default ax_ratio=50 is more than enough for producing consistent results.\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.discretize-Union{Tuple{ftype}, Tuple{PlaneFaultDomain{ftype,2,T} where T,Any,Any}} where ftype<:JuEQ.PlaneFault",
    "page": "Public",
    "title": "JuEQ.discretize",
    "category": "method",
    "text": "discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ; buffer=:auto) where {ftype <: PlaneFault}\n\nGenerate the grid for given 2D fault domain.\n\nArguments\n\nΔx, Δξ: grid space along-strike and along-downdip respectively\nbuffer::Union{T, Symbol}: length of buffer size for introducing zero-dislocation area at along-strike edges of defined fault domain.\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.fault-Union{Tuple{T}, Tuple{N}, Tuple{Type{#s32} where #s32<:PlaneFault,T,Tuple{Vararg{T,N}}}} where T where N",
    "page": "Public",
    "title": "JuEQ.fault",
    "category": "method",
    "text": "fault(ftype::Type{<:PlaneFault}, dip, span)\n\nGenerate a fault given the fault type, dip angle and its spatial span.\n\nArguments\n\nftype::Type{<:PlaneFault}: type of plane fault\ndip: dip angle in degree\nspan: spatial span of fault size\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.friction-Union{Tuple{T}, Tuple{CForm,T,T,T,T,T,T,T}} where T<:Number",
    "page": "Public",
    "title": "JuEQ.friction",
    "category": "method",
    "text": "friction(::FrictionLawForm, v::T, θ::T, L::T, a::T, b::T, f0::T, v0::T) where {T<:Number}\n\nCalculate friction given by the form of fomula as well as other necessary parameters.\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.properties-Union{Tuple{dim}, Tuple{PlaneFaultDomain,BoundaryElementGrid{dim,isuniform,T} where T where isuniform}} where dim",
    "page": "Public",
    "title": "JuEQ.properties",
    "category": "method",
    "text": "properties(fa::PlaneFaultDomain, gd::BoundaryElementGrid{dim}; _kwargs...) where {dim}\n\nEstablishing a material-properties-profile given by the fault domain and grids. User must provide the     necessary parameters in according to the grid size specified or just a scalar for broadcasting.\n\nArguments that are needed:\n\na: contrib from velocity.\nb: contrib from state.\nL: critical distance.\nσ: effective normal stress.\nvpl: plate rate.\nf0: ref. frictional coeff.\nv0: ref. velocity.\n\nArguments that are optional\n\nk: stiffness tensor. If :auto, it will automatically calculate by seeking λ and μ otherwise should be a valid file path to a *.jld2.\nη: radiation damping. If :auto, it will automatically seek μ and vs and use μ  2vs.\nvs: shear wave velocity.\nλ: Lamé\'s first parameter\nμ: shear modulus\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#JuEQ.stiffness_tensor-Tuple{PlaneFaultDomain,JuEQ.BoundaryElementGrid,JuEQ.HomogeneousElasticProperties}",
    "page": "Public",
    "title": "JuEQ.stiffness_tensor",
    "category": "method",
    "text": "stiffness_tensor(fa::PlaneFaultDomain, gd::BoundaryElementGrid, ep::HomogeneousElasticProperties)\n\nCalculate the reduced stiffness tensor. For 2D fault, the final result will be dimensionally reduced to a 3D array     due to the translational & reflective & perodic symmetry, such that the tensor contraction will be equivalent to convolution,     hence we could use FFT for better performace.\n\nNote\n\nFaults are originated from surface and extends downwards, thus dep = 0\n\n\n\n\n\n"
},

{
    "location": "public_interface.html#Interfaces-1",
    "page": "Public",
    "title": "Interfaces",
    "category": "section",
    "text": "Modules = [JuEQ]\nPrivate = false\nOrder = [:type, :function, :constant, :macro]"
},

{
    "location": "private_interface.html#",
    "page": "Private",
    "title": "Private",
    "category": "page",
    "text": ""
},

{
    "location": "private_interface.html#Private-Interface-1",
    "page": "Private",
    "title": "Private Interface",
    "category": "section",
    "text": ""
},

{
    "location": "private_interface.html#Index-1",
    "page": "Private",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"private_interface.md\"]"
},

{
    "location": "private_interface.html#JuEQ.HomogeneousElasticProperties",
    "page": "Private",
    "title": "JuEQ.HomogeneousElasticProperties",
    "category": "type",
    "text": "Okada\'s dc3d only applies on isotropic materials,     therefore, elastic modulus are constrained to be scalars.\n\n\n\n\n\n"
},

{
    "location": "private_interface.html#JuEQ.TmpVariable",
    "page": "Private",
    "title": "JuEQ.TmpVariable",
    "category": "type",
    "text": "Temporal variable in solving ODEs aimed to avoid allocation overheads.\n\n\n\n\n\n"
},

{
    "location": "private_interface.html#JuEQ.applied_unit_dislocation-Tuple{Type{NormalFault}}",
    "page": "Private",
    "title": "JuEQ.applied_unit_dislocation",
    "category": "method",
    "text": "For noraml fault, it should of course be [0., -1., 0.]. However, in term of force balance, it is quivalent to thrust fault if dip angle are constrained within [0, π/2] in fact.\n\nThe unit of unit dislocation below is the same of v * t at set by user so to avoid normalization step.\n\n\n\n\n\n"
},

{
    "location": "private_interface.html#JuEQ.shear_traction-Tuple{Type{#s32} where #s32<:JuEQ.AbstractFault}",
    "page": "Private",
    "title": "JuEQ.shear_traction",
    "category": "method",
    "text": "shear_traction(::Type{<:PlaneFault}, u, λ, μ, dip)\n\nCalculate the shear traction on the fault plane w.r.t. fault types.\n\nArguments\n\nu::AbstractArray{<:Number, 1}: the output from dc3d_okada\nλ::Number: Lamé\'s first parameter\nμ::Number: shear modulus\ndip::Number: plane dip angle\n\nReference\n\nA good reference is at Displacement & Strain & Stress.\n\n\n\n\n\n"
},

{
    "location": "private_interface.html#JuEQ.stiffness_periodic_boundary_condition-Union{Tuple{T}, Tuple{T,T,T,T,T,T,AbstractArray{T,1},AbstractArray{T,1},AbstractArray{T,1},Integer,T}} where T<:Number",
    "page": "Private",
    "title": "JuEQ.stiffness_periodic_boundary_condition",
    "category": "method",
    "text": "Periodic boundary condition for 2D faults.\n\nArguments\n\nsame as dc3d_okada, see dc3d for details.\nax::T: along-strike fault length\nnrept::Integer: (half) number of repetition, as denoted by -npret: nrept\nlrept::T: length of repetition interval, see Note below\n\nNote\n\nThe buffer block is evenly distributed on the two along-strike edges, each of which contains half of that.\n\n\n\n\n\n"
},

{
    "location": "private_interface.html#Interfaces-1",
    "page": "Private",
    "title": "Interfaces",
    "category": "section",
    "text": "Modules = [JuEQ]\nPublic = false\nOrder = [:type, :function, :constant, :macro]"
},

]}
