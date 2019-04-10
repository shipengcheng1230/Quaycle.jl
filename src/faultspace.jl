## fault space

export fault

abstract type FaultSpace end

struct CentralSymmetryFS{G<:SimpleMesh, FT<:PlaneFault, T<:Real} <: FaultSpace
    mesh::G
    dip::T
    faulttype::FT
end

function fault(::Val{:CSFS}, ft::FT, dip::T, span::T, dξ::T) where {T<:Real, FT<:PlaneFault}
    geo = SimpleLine((0., 0., 0.,), (span,))
    mesh = gen_mesh(Val(:SimpleGrid), geo, dξ)
    return CentralSymmetryFS(mesh, dip, ft)
end

function fault(::Val{:CSFS}, ft::FT, dip::T, span::NTuple{2,T}, dxdξ::NTuple{2,T}) where {T<:Real, FT<:PlaneFault}
    geo = SimpleRect((0., 0., 0.,), span)
    mesh = gen_mesh(Val(:SimpleGrid), geo, dxdξ...)
    return CentralSymmetryFS(mesh, dip, ft)
end


fault(fst::Val{:CSFS}, ft::FT, dip, span, dξ) where {FT<:PlaneFault} = fault(fst, ft, dip, promote(span, dξ)...)
fault(fst::Val{:CSFS}, ft::FT, dip, span1, span2, dx, dξ) where {FT<:PlaneFault} = fault(fst, ft, promote(dip, span1, span2, dx, dξ)...)
fault(fst::Val{:CSFS}, ft::FT, dip::T, span1::T, span2::T, dx::T, dξ::T) where {T<:Real, FT<:PlaneFault} = fault(fst, ft, dip, (span1, span2), (dx, dξ))

fault(fst::Val{:CSFS}, ft::STRIKING, span, dξ) = fault(fst, ft, 90.0, span, dξ)
fault(fst::Val{:CSFS}, ft::STRIKING, span1, span2, dx, dξ) = fault(fst, ft, 90.0, span1, span2, dx, dξ)

ft1 = fault(Val(:CSFS), STRIKING(), 90.0, 10.0, 1.0)
ft2 = fault(Val(:CSFS), STRIKING(), 10, 1.0)

ft3 = fault(Val(:CSFS), STRIKING(), 90.0, (10.0, 5.0), (1.0, 0.5))
ft4 = fault(Val(:CSFS), STRIKING(), 90.0, 10.0, 5.0, 1.0, 0.5)
ft5 = fault(Val(:CSFS), STRIKING(), 90.0, 10, 5, 1, 0.5)
ft6 = fault(Val(:CSFS), STRIKING(), 10.0, 5.0, 1.0, 0.5)

ft1 == ft2

ft1.mesh.aξ == ft2.mesh.aξ
