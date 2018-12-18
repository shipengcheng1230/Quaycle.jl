var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#JuEQ.jl-Documentation-1",
    "page": "Home",
    "title": "JuEQ.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "This is a suite for numerically simulating earthquake sequences in Julia. The purpose of this package is to provide efficient Julia implementations for simulations in the field of earthquake physics.Features of this package are listed as below:Rate-State Friction Law\nOkada\'s Dislocation Method\nBoundary Element Method (Quasi-dynamic)Features to be implemented:Viscoelastic Relaxation (priority)\nFully Elastodynamic Effect\nOff-Fault Materials Effect"
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Get the latest version with Julia\'s package manager:] add https://github.com/shipengcheng1230/JuEQ.jlTo load the package:using JuEQ"
},

{
    "location": "#Acknowledgements-1",
    "page": "Home",
    "title": "Acknowledgements",
    "category": "section",
    "text": "The simulation of episodic seismic and slow slip events using boundary-element-method is largely benifit from Yajing Liu original Fortran code."
},

{
    "location": "quasi_dynamic_intro/#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "quasi_dynamic_intro/#Quasi-dynamic-Simulation-1",
    "page": "Introduction",
    "title": "Quasi-dynamic Simulation",
    "category": "section",
    "text": ""
},

{
    "location": "quasi_dynamic_intro/#Basic-Theory-1",
    "page": "Introduction",
    "title": "Basic Theory",
    "category": "section",
    "text": "The governing equation is that at every time step, shear stress across the fault plane equals to frictional force plus a radiation damping term for approximating wave propagation effect:τ = σf + ηVHere μ is shear stress across the fault plain. Using Okada\'s dislocation theory, it can be shown as:τ = K  δwhere K is the so-called stiffness tensor, depicting relationship between displacements at one position regarding to dislocations somewhere else. δ is the dislocation, i.e. displacement at everywhere on the fault.  denotes tensor contraction.Back to f, we use rate-and-state frictional law to calculate its value, specifically as below:f(V θ) = f_0 + a lnfracVV_0 + b lnleft(fracV_0 θLright)where f_0 and V_0 are reference friction coefficient and velocity, V and θ are velocity and state variable based on which frictional force is. a and b are two frictional parameters denoting contributions each of which comes from velocity and state variable respectively. L is critical distance after which frictional force return to new steady state.Sometimes people use regularized form to avoid infinity when V  0, namely:f(V θ) = a sinh ^-1leftfracV2V_0 expfracf_0 + b lnleft(V_0 θLright)arightThere are many state evolution law that describes how state variable θ changes with time, one of which that most widely used is Dieterich law:fracmathrmdθmathrmdt = 1 - fracV θLFurther, η is a damping coefficient whose value is often chosen as μ  2mathrmVs where μ is shear modulus and mathrmVs shear wave velocity and σ is the effective normal stress.To simulate how fault evolves with time, we then take the derivative of the governing equation:fracmathrmd τmathrmd t = fracmathrmd f(V θ)mathrmd t + η fracmathrmd Vmathrmd t\n= fracmathrmd fmathrmd V fracmathrmd Vmathrmd t + fracmathrmd fmathrmd θ fracmathrmd θmathrmd t + η fracmathrmd Vmathrmd tThus we arrive at:fracmathrmd Vmathrmd t = fracfracmathrmd τmathrmd t - fracmathrmd fmathrmd θ fracmathrmd θmathrmd tfracmathrmd fmathrmd V + ηwhere fracmathrmd τmathrmd t = K  (mathrmV_pl - V) where mathrmV_pl is the plate rate.note: Note\nThe direction of relative velocity, namely mathrmV_pl - V, must be in accordance to the direction of K which, here, we use the same meaning as Rice, J. (1993).Hence, with both derivatives of velocity V and state variable θ, we are able to discover how fault evolves with various parameters settings."
},

{
    "location": "generated/bp1/#",
    "page": "Example 1: 1D fault",
    "title": "Example 1: 1D fault",
    "category": "page",
    "text": "EditURL = \"https://github.com/shipengcheng1230/JuEQ.jl/blob/master/examples/bp1.jl\"note: Note\nThis example is from Benchmark Problem 1 (hence referred as BP1)."
},

{
    "location": "generated/bp1/#Define-parameters-1",
    "page": "Example 1: 1D fault",
    "title": "Define parameters",
    "category": "section",
    "text": "First, we load the packageusing JuEQ\nusing PlotsInstead of using SI unit, we refactor ours into the follow:ms2mmyr = 365 * 86400 * 1e3 # convert velocity from m/s to mm/yr\nρ = 2670.0 # density [kg/m³]\nvs = 3464.0 # shear wave velocity [m/s]\nσ = 500.0 # effective normal stress [bar]\na0 = 0.010 # frictional paramter `a` in vw zone\namax = 0.025 # frictional paramter `a` in vs zone\nb0 = 0.015 # frictional paramter `b`\nL = 8.0 # critical distance [mm]\nvpl = 1e-9 * ms2mmyr # plate rate [mm/yr]\nvinit = 1e-9 * ms2mmyr # initial velocity [mm/yr]\nv0 = 1e-6 * ms2mmyr # reference velocity [mm/yr]\nf0 = 0.6 # reference frictional coefficient\nH = 15.0 # vw zone [km]\nh = 3.0 # vw-vs changing zone [km]\nWf = 40.0 # fault depth [km]\nΔz = 100.0e-3 # grid size interval [km]\ntf = 400.0; # simulation time [yr]warning: Warning\nMake sure your units are consistent across the whole variable space. Pontenial imporvement may incoporate Unitful.jl package.Then we arrive at some parameters that are implicit by above:μ = vs^2 * ρ / 1e5 / 1e6 # shear modulus [bar·km/mm]\nλ = μ # poisson material\nη = μ / 2(vs * 1e-3 * 365 * 86400)\nngrid = round(Int, Wf / Δz); # number of gridsNow, we start to construct our model using parameters above. First, we create a \'fault\' by specifying fault type and depth:tip: Tip\nHere, we do not need to provide dip for strike-slip fault as it automatically choose 90. See fault."
},

{
    "location": "generated/bp1/#Construct-Model-1",
    "page": "Example 1: 1D fault",
    "title": "Construct Model",
    "category": "section",
    "text": "fa = fault(StrikeSlipFault, Wf);Next, we generate the grid regarding the fault we just created by giving number of grids:note: Note\nThis package use ξ for denoting downdip coordinate and x for along-strike one. See discretize.gd = discretize(fa; nξ=ngrid);Next, we construct the required frictional parameter profile:z = -gd.ξ\naz = fill(a0, size(z))\naz[z .≥ (H + h)] .= amax\naz[H .< z .< H + h] = a0 .+ (amax - a0) / (h / Δz) * collect(1: Int(h / Δz));Then, we provide the required initial condition satisfying uniform slip distribution over the depth:τ0 = σ * amax * asinh(vinit / 2v0 * exp((f0 + b0 * log(v0 / vinit)) / amax)) + η * vinit\nτz = fill(τ0, size(z))\nθz = @. L / v0 * exp(az / b0 * log(2v0 / vinit * sinh((τz - η * vinit) / az / σ)) - f0 / b0)\nvz = fill(vinit, size(z))\nu0 = hcat(vz, θz);Let\'s simulate only the first 200 years:tspan = (0., 200.);Finally, we provide the material properties w.r.t. our \'fault\', \'grid\' as well as other necessary parameters predefined using the same grid size & dimension:mp = properties(;fault=fa, grid=gd, parameters=[:a=>az, :b=>b0, :L=>L, :σ=>σ, :η=>η, :k=>[:λ=>λ, :μ=>μ], :vpl=>vpl, :f0=>f0, :v0=>v0]);tip: Tip\nCheck properties for extended options.Check our profile now:plot([mp.a, mp.b], z, label=[\"a\", \"b\"], yflip=true, ylabel=\"Depth (km)\")We then contruct the ODEProblem as following by stating which state evolution law to use and frcitonal law form, plus initial condition and simulation time:prob = EarthquakeCycleProblem(gd, mp, u0, tspan; se=DieterichStateLaw(), fform=RForm());"
},

{
    "location": "generated/bp1/#Solve-Model-1",
    "page": "Example 1: 1D fault",
    "title": "Solve Model",
    "category": "section",
    "text": "We then solve the ODEs:sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6);tip: Tip\nFor details of solving options, see here.tip: Tip\nRaise the accuracy option if you get instability when solving these ODEs."
},

{
    "location": "generated/bp1/#Results-1",
    "page": "Example 1: 1D fault",
    "title": "Results",
    "category": "section",
    "text": "The first event happens at around 196 year:maxv = max_velocity(sol)\nplot(sol.t, log10.(maxv / ms2mmyr), xlabel=\"Time (year)\", ylabel=\"Max Velocity (log10 (m/s))\", xlims=(190, 200), label=\"\")note: Note\nClick here for the slip evolution over 3000 years simulation. It may need some time to load the page.This page was generated using Literate.jl."
},

{
    "location": "generated/otfsync/#",
    "page": "Example 2: 2D fault",
    "title": "Example 2: 2D fault",
    "category": "page",
    "text": "EditURL = \"https://github.com/shipengcheng1230/JuEQ.jl/blob/master/examples/otfsync.jl\"note: Note\nThis example is adapted from Wei, 2016 AGUtip: Tip\nIt will automatically use parallel scheme if nprocs() ≂̸ 1 when building stiffness tensor. To do so:using Distributed\naddprocs(4); # add # of cores you desire\nusing JuEQFirst, we load the package and define some basic parameters:using JuEQ\nusing Plots\n\nms2mmyr = 365 * 86400 * 1e3\nρ = 2670.0 # kg/m³\ncs = 3044.0 # m/s\nvpl = 100.0 # mm/yr\nv0 = 3.2e4 # mm/yr\nf0 = 0.6;Then we come to parameters implicit by above:μ = 0.3 # Bar·km/mm\nλ = μ # poisson material\nα = (λ + μ) / (λ + 2μ)\nη = μ / 2(cs * 1e-3 * 365 * 86400); # Bar·yr/mmCreate a fault:fa = fault(StrikeSlipFault, (80., 10.));Generate grids:gd = discretize(fa; nx=160, nξ=20, buffer_ratio=1);tip: Tip\nIt is recommended (from Yajing Liu\'s personal communication) to add buffer zones adjacent the horizontal edges to immitate zero dislocation at the ridge region. Basically, it affects how the stiffness tensor are periodically summed. To what extent it alters the results remains further testing.Under the hood, it shall impose buffer areas on both sides of along-strike, each of which has a length of bufferratio/2*fa[:x]. Thus, the stiffness contributions falling into those buffer zone shall be neglected, which is equivalent to impose zero-slip correspondingly.Time for us to establish frictional parameters profile:a = 0.015 .* ones(gd.nx, gd.nξ)\nb = 0.0115 .* ones(gd.nx, gd.nξ)\nleft_patch = @. -25. ≤ gd.x ≤ -5.\nright_patch = @. 5. ≤ gd.x ≤ 25.\nvert_patch = @. -6. ≤ gd.z ≤ -1.\nb[xor.(left_patch, right_patch), vert_patch] .= 0.0185\namb = a - b\nσmax = 500.\nσ = [min(σmax, 15. + 180. * z) for z in -gd.z]\nσ = Matrix(repeat(σ, 1, gd.nx)\')\nL = 12.;Check our profile:p1 = heatmap(amb\',\n    xticks=(0: 10/gd.Δx: gd.nx, -fa.span[1]/2: 10: fa.span[1]/2),\n    yticks=(0: 5/gd.Δξ: gd.nξ, 0: -5: -fa.span[2]),\n    yflip=true, color=:isolum, aspect_ratio=2, title=\"a-b\"\n    );\n\np2 = heatmap(σ\',\n    xticks=(0: 10/gd.Δx: gd.nx, -fa.span[1]/2: 10: fa.span[1]/2),\n    yticks=(0: 5/gd.Δξ: gd.nξ, 0: -5: -fa.span[2]),\n    yflip=true, color=:isolum, aspect_ratio=2, title=\"\\\\sigma\"\n    );\n\nplot(p1, p2, layout=(2, 1))Construct our material property profile:mp = properties(fa, gd, [:a=>a, :b=>b, :L=>L, :σ=>σ, :η=>η, :k=>[:λ=>λ, :μ=>μ], :vpl=>vpl, :f0=>f0, :v0=>v0]);Provide the initial condition:vinit = vpl .* ones(gd.nx, gd.nξ)\nθ0 = L ./ vinit ./ 1.1\nu0 = cat(vinit, θ0, dims=3);Get our ODEs problem:prob = EarthquakeCycleProblem(gd, mp, u0, (0., 18.); se=DieterichStateLaw(), fform=CForm());Solve the model:sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6);Take a look at the max velocity:maxv = max_velocity(sol)\nplot(sol.t, log10.(maxv / ms2mmyr), xlabel=\"Time (year)\", ylabel=\"Max Velocity (log10 (m/s))\", label=\"\")View some snapshots to see the rupture (quasi-dynamic) patterns:ind = argmax(maxv)\nmyplot = (ind) -> heatmap(log10.(sol.u[ind][:,:,1]./ms2mmyr)\',\n    xticks=(0: 10/gd.Δx: gd.nx, -fa.span[1]/2: 10: fa.span[1]/2),\n    yticks=(0: 5/gd.Δξ: gd.nξ, 0: -5: -fa.span[2]),\n    yflip=true, color=:isolum, aspect_ratio=2, title=\"t = $(sol.t[ind])\")\n\nsnaps = [myplot(i) for i in ind-700: 200: ind+500]\n\nplot(snaps..., layout=(length(snaps), 1), size=(600, 1800))This page was generated using Literate.jl."
},

{
    "location": "public_interface/#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "public_interface/#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": ""
},

{
    "location": "public_interface/#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public_interface.md\"]"
},

{
    "location": "public_interface/#JuEQ.DieterichStateLaw",
    "page": "Public",
    "title": "JuEQ.DieterichStateLaw",
    "category": "type",
    "text": "fracmathrmdθmathrmdt = 1 - fracV θL\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.PrzStateLaw",
    "page": "Public",
    "title": "JuEQ.PrzStateLaw",
    "category": "type",
    "text": "fracmathrmdθmathrmdt = 1 - (fracV θ2L)^2\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.RuinaStateLaw",
    "page": "Public",
    "title": "JuEQ.RuinaStateLaw",
    "category": "type",
    "text": "fracmathrmdθmathrmdt = -fracV θL * logfracV θL\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.DECallbackSaveToFile-Tuple{IOStream,IOStream}",
    "page": "Public",
    "title": "JuEQ.DECallbackSaveToFile",
    "category": "method",
    "text": "DECallbackSaveToFile(iot::IOStream, iou::IOStream)\n\nConstruct a functional callback to write ODESolution (t & u) into file. The reason to separate t and u is for more easily reshape u w.r.t grids specification. It right now falls on users\' memory on what the type of solution is for accurately retrieving results.\n\nArguments\n\niot::IOStream: stream pointing to solution of time\niou::IOStream: stream pointing to solution of domain\n\nNote It is strongly not recommended to use \"skipping\" scheme (by defining thrd and dts(a) for each case) when solution is too oscillated.\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.EarthquakeCycleProblem-Union{Tuple{dim}, Tuple{BoundaryElementGrid,PlaneMaterialProperties,AbstractArray,Tuple{Vararg{T,N}} where T where N}} where dim",
    "page": "Public",
    "title": "JuEQ.EarthquakeCycleProblem",
    "category": "method",
    "text": "EarthquakeCycleProblem(p::PlaneMaterialProperties, u0, tspan; se=DieterichStateLaw(), fform=CForm())\n\nReturn an ODEProblem that encapsulate all the parameters and functions required for simulation. For the entailing usage, please refer DifferentialEquations.jl\n\nArguments\n\ngd::BoundaryElementGrid: grids for fault domain.\np::PlaneMaterialProperties: material profile.\nu0::AbstractArray: initial condition, should be organized such that the first of last dim is velocity while the 2nd of last dim is state.\ntspan::NTuple: time interval to be simulated.\nse::StateEvolutionLaw: state evolution law to be applied.\nfform::FrictionLawForm: forms of frictional law to be applied.\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.dc3d_okada-Union{Tuple{A}, Tuple{T}, Tuple{T,T,T,T,T,T,Union{SubArray, A},Union{SubArray, A},A}} where A<:Union{AbstractArray{T,1}, AbstractArray{T,2}} where T<:Number",
    "page": "Public",
    "title": "JuEQ.dc3d_okada",
    "category": "method",
    "text": "Calculate displacements and gradient of displacements due to a dislocation in an elastic isotropic halfspace. See dc3d for details.\n\ntest/test_okada.dat is obtained using DC3dfortran\n\nAn example wrapper for DC3D in julia as below:\n\nfunction dc3d_fortran(x::T, y::T, z::T, α::T, dep::T, dip::T, al1::T, al2::T, aw1::T, aw2::T,\n    disl1::T, disl2::T, disl3::T) where {T <: AbstractFloat}\n\n    # initial return values\n    # `RefValue{T}` may be also viable other than `Array{T, 1}`\n    ux = Array{Float64}(1)\n    uy = Array{Float64}(1)\n    uz = Array{Float64}(1)\n    uxx = Array{Float64}(1)\n    uyx = Array{Float64}(1)\n    uzx = Array{Float64}(1)\n    uxy = Array{Float64}(1)\n    uyy = Array{Float64}(1)\n    uzy = Array{Float64}(1)\n    uxz = Array{Float64}(1)\n    uyz = Array{Float64}(1)\n    uzz = Array{Float64}(1)\n    iret = Array{Int64}(1)\n\n    # call okada\'s code which is renamed as \"__dc3d__\" (see binding rename shown below)\n    # input args tuple must be syntactically written instead of a variable assigned\n    # macros could be used to simplify this in the future\n    ccall((:__dc3d__, \"dc3d.so\"), Void,\n        (\n            Ref{Float64},\n            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},\n            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},\n            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},\n            Ref{Int64},\n        ),\n        α, x, y, z, dep, dip, al1, al2, aw1, aw2, disl1, disl2, disl3,\n        ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz,\n        iret,\n    )\n\n    # results valid iff iret[1] == 0\n    return (\n        iret[1],\n        ux[1], uy[1], uz[1],\n        uxx[1], uyx[1], uzx[1],\n        uxy[1], uyy[1], uzy[1],\n        uxz[1], uyz[1], uzz[1]\n    )\nend\n\nThe corresponding fortran module is:\n\nMODULE okada\n  USE, INTRINSIC :: iso_c_binding\n  IMPLICIT NONE\nCONTAINS\n\n  SUBROUTINE dc3d_wrapper(&\n       & alpha, &\n       & x, y, z, &\n       & depth, dip, &\n       & al1, al2, &\n       & aw1, aw2, &\n       & disl1, disl2, disl3, &\n       & ux, uy, uz, &\n       & uxx, uyx, uzx, &\n       & uxy, uyy, uzy, &\n       & uxz, uyz, uzz, &\n       & iret) BIND(C, NAME=\'__dc3d__\')\n\n    REAL*8 :: &\n         & alpha, &\n         & x, y, z, &\n         & depth, dip, &\n         & al1, al2, &\n         & aw1, aw2, &\n         & disl1, disl2, disl3, &\n         & ux, uy, uz, &\n         & uxx, uyx, uzx, &\n         & uxy, uyy, uzy, &\n         & uxz, uyz, uzz\n\n    INTEGER*8 :: iret\n\n    CALL dc3d(&\n         & alpha, &\n         & x, y, z, &\n         & depth, dip, &\n         & al1, al2, &\n         & aw1, aw2, &\n         & disl1, disl2, disl3, &\n         & ux, uy, uz, &\n         & uxx, uyx, uzx, &\n         & uxy, uyy, uzy, &\n         & uxz, uyz, uzz, &\n         & iret)\n\n  END SUBROUTINE dc3d_wrapper\n\nEND MODULE okada\n\nA sample of makefile is as below:\n\n# Build Okada\'s code for calculating deformation due to a fault model\n#\nCC = gfortran\nCFLAGS = -fPIC -w -O3\nLDFLAGS = -shared\n\nSRCS = dc3d.f okada.f90\nOBJS = $(SRCS:.c=.o)\n\nTARGET = dc3d.so\n\n$(TARGET): $(OBJS)\n    $(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS)\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.discretize-Union{Tuple{T}, Tuple{ftype}, Tuple{PlaneFaultDomain{ftype,1,T},T}} where T<:Number where ftype<:JuEQ.PlaneFault",
    "page": "Public",
    "title": "JuEQ.discretize",
    "category": "method",
    "text": "discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ::T; ax_ratio=12.5)\n\nGenerate the grid for given 1D fault domain. The grids will be forced to start at (x=0, y=0, z=0).\n\nArguments\n\nΔξ: grid space along-downdip\nax_ratio::Number: ration of along-strike length agsinst along-downdip length for mimicing an extended   2d (x & ξ) fault represented by 1d (ξ) domain. Default ax_ratio=12.5 is more than enough for producing consistent results.\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.discretize-Union{Tuple{T}, Tuple{ftype}, Tuple{PlaneFaultDomain{ftype,2,T},T,T}} where T<:Number where ftype<:JuEQ.PlaneFault",
    "page": "Public",
    "title": "JuEQ.discretize",
    "category": "method",
    "text": "discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ; buffer=:auto) where {ftype <: PlaneFault}\n\nGenerate the grid for given 2D fault domain. The grids will be forced to start at (z=0) and spread symmetrically along x-axis w.r.t y-z plane.     By such setting, we would be able to utilize the symmetry properties of stiffness tensor for performance speed up.\n\nArguments\n\nΔx, Δξ: grid space along-strike and along-downdip respectively\n`buffer_ratio::Number: ration of buffer size against along-strike length for introducing zero-dislocation area at along-strike edges of defined fault domain.\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.fault-Union{Tuple{T}, Tuple{N}, Tuple{Type{#s33} where #s33<:PlaneFault,T,Tuple{Vararg{T,N}}}} where T where N",
    "page": "Public",
    "title": "JuEQ.fault",
    "category": "method",
    "text": "fault(ftype::Type{<:PlaneFault}, dip, span)\n\nGenerate a fault given the fault type, dip angle and its spatial span.\n\nArguments\n\nftype::Type{<:PlaneFault}: type of plane fault\ndip: dip angle in degree\nspan: spatial span of fault size\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.friction-Union{Tuple{T}, Tuple{CForm,T,T,T,T,T,T,T}} where T<:Number",
    "page": "Public",
    "title": "JuEQ.friction",
    "category": "method",
    "text": "friction(::FrictionLawForm, v::T, θ::T, L::T, a::T, b::T, f0::T, v0::T) where {T<:Number}\n\nCalculate friction given by the form of fomula as well as other necessary parameters.\n\nConventional Form:\n\nf(V θ) = f_0 + a lnfracVV_0 + b lnleft(fracV_0 θLright)\n\nRegularized Form:\n\nf(V θ) = a sinh ^-1leftfracV2V_0 expfracf_0 + b lnleft(V_0 θLright)aright\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.max_velocity-Tuple{AbstractArray{T,1} where T,AbstractArray,Function}",
    "page": "Public",
    "title": "JuEQ.max_velocity",
    "category": "method",
    "text": "max_velocity(t::AbstractVector, u::AbstractArray, getu::Function)\n\nReturn max velocity across the fault at each time step. A number of convenient interfaces for common output are implemented.\n\nArguments\n\nt::AbstractVector: vector of time steps\nu::AbstractArray: array of solution\ngetu::Function: method for retrieving velocity section at each time step\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.moment_magnitude-Union{Tuple{T}, Tuple{T,T,T}} where T<:Number",
    "page": "Public",
    "title": "JuEQ.moment_magnitude",
    "category": "method",
    "text": "Calculate moment magnitude.\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.properties-Union{Tuple{dim}, Tuple{ftype}, Tuple{PlaneFaultDomain{ftype,dim,T} where T,BoundaryElementGrid,AbstractArray{#s46,N} where N where #s46<:Pair}} where dim where ftype<:JuEQ.PlaneFault",
    "page": "Public",
    "title": "JuEQ.properties",
    "category": "method",
    "text": "properties(fa::PlaneFaultDomain, gd::BoundaryElementGrid{dim}; _kwargs...) where {dim}\n\nEstablishing a material-properties-profile given by the fault domain and grids. User must provide the     necessary parameters in according to the grid size specified or just a scalar for broadcasting.\n\nArguments that are required:\n\na: contrib from velocity.\nb: contrib from state.\nL: critical distance.\nσ: effective normal stress.\nη: radiation damping. It is recommended to set as μ  2mathrmVs where μ is shear modulus and mathrmVs shear wave velocity.\nvpl: plate rate.\nf0: ref. frictional coeff.\nv0: ref. velocity.\n\nArguments that need options\n\nk: stiffness tensor.\n(1) Providing shear modulus denoted as μ and Lamé\'s first parameter denoted as λ (same as μ if missing),  then calculate it based on grid and fault domain, choosing parallel scheme if nprocs() != 1.  (2) A valid file path to a .jld2 that contains valid stiffness tensor. No verification will be performed here.  (3) an AbstractArray represent the pre-calculated stiffness tensor. No verification will be performed here.\n\n\n\n\n\n"
},

{
    "location": "public_interface/#JuEQ.stiffness_tensor-Union{Tuple{T}, Tuple{ftype}, Tuple{PlaneFaultDomain{ftype,1,T},BoundaryElementGrid{1},HomogeneousElasticProperties}} where T<:Number where ftype<:JuEQ.PlaneFault",
    "page": "Public",
    "title": "JuEQ.stiffness_tensor",
    "category": "method",
    "text": "stiffness_tensor(fa::PlaneFaultDomain, gd::BoundaryElementGrid, ep::HomogeneousElasticProperties)\n\nCalculate the reduced stiffness tensor. For 2D fault, the final result will be dimensionally reduced to a 3D array     due to the translational & reflective & perodic symmetry, such that the tensor contraction will be equivalent to convolution,     hence we could use FFT for better performace.\n\nNote\n\nFaults are originated from surface and extends downwards, thus dep = 0\n\n\n\n\n\n"
},

{
    "location": "public_interface/#Interfaces-1",
    "page": "Public",
    "title": "Interfaces",
    "category": "section",
    "text": "Modules = [JuEQ]\nPrivate = false\nOrder = [:type, :function, :constant, :macro]"
},

{
    "location": "private_interface/#",
    "page": "Private",
    "title": "Private",
    "category": "page",
    "text": ""
},

{
    "location": "private_interface/#Private-Interface-1",
    "page": "Private",
    "title": "Private Interface",
    "category": "section",
    "text": ""
},

{
    "location": "private_interface/#Index-1",
    "page": "Private",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"private_interface.md\"]"
},

{
    "location": "private_interface/#JuEQ.HomogeneousElasticProperties",
    "page": "Private",
    "title": "JuEQ.HomogeneousElasticProperties",
    "category": "type",
    "text": "Okada\'s dc3d only applies on isotropic materials,     therefore, elastic modulus are constrained to be scalars.\n\n\n\n\n\n"
},

{
    "location": "private_interface/#JuEQ.ODEStateVariable",
    "page": "Private",
    "title": "JuEQ.ODEStateVariable",
    "category": "type",
    "text": "Intermediate variable in solving ODEs aimed to avoid allocation overheads.\n\n\n\n\n\n"
},

{
    "location": "private_interface/#JuEQ.applied_unit_dislocation-Tuple{Type{NormalFault}}",
    "page": "Private",
    "title": "JuEQ.applied_unit_dislocation",
    "category": "method",
    "text": "For noraml fault, it should of course be [0., -1., 0.]. However, in term of force balance, it is quivalent to thrust fault if dip angle are constrained within [0, π/2] in fact.\n\nThe unit of unit dislocation below is the same of v * t at set by user so to avoid normalization step.\n\n\n\n\n\n"
},

{
    "location": "private_interface/#JuEQ.shear_traction-Tuple{Type{#s46} where #s46<:JuEQ.AbstractFault}",
    "page": "Private",
    "title": "JuEQ.shear_traction",
    "category": "method",
    "text": "shear_traction(::Type{<:PlaneFault}, u, λ, μ, dip)\n\nCalculate the shear traction on the fault plane w.r.t. fault types.\n\nArguments\n\nu::AbstractArray{<:Number, 1}: the output from dc3d_okada\nλ::Number: Lamé\'s first parameter\nμ::Number: shear modulus\ndip::Number: plane dip angle\n\nReference\n\nA good reference is at Displacement & Strain & Stress.\n\n\n\n\n\n"
},

{
    "location": "private_interface/#JuEQ.stiffness_periodic_boundary_condition!-Union{Tuple{T}, Tuple{AbstractArray{T,1},T,T,T,T,T,T,AbstractArray{T,1},AbstractArray{T,1},AbstractArray{T,1},Integer,T}} where T<:Number",
    "page": "Private",
    "title": "JuEQ.stiffness_periodic_boundary_condition!",
    "category": "method",
    "text": "Periodic boundary condition for 2D faults.\n\nArguments\n\nu::AbstractVector: In-place output which is a 12-elements vector (exactly the output of dc3d_okada). No assertion here imposed.\nsame as dc3d_okada, see dc3d for details.\nnrept::Integer: (half) number of repetition, as denoted by -npret: nrept\nlrept::Number: length of repetition interval, see Note below\n\nNote\n\nThe buffer block length is (buffer_ratio - 1) multipled by along-strike length.\n\n\n\n\n\n"
},

{
    "location": "private_interface/#Interfaces-1",
    "page": "Private",
    "title": "Interfaces",
    "category": "section",
    "text": "Modules = [JuEQ]\nPublic = false\nOrder = [:type, :function, :constant, :macro]"
},

]}
