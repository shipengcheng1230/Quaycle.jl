export vtu_output

const gmshcelltype2vtkcelltype = Dict(
    1 => VTKCellTypes.VTK_LINE,
    3 => VTKCellTypes.VTK_QUAD,
    4 => VTKCellTypes.VTK_TETRA,
    5 => VTKCellTypes.VTK_HEXAHEDRON,
)

abstract type ParaviewOutputCache end

struct VTUStructuredCellDataCache{D, V, A, C, P, L} <: ParaviewOutputCache
    tagmap::D
    etag::V
    dat::A
    cells::C
    pts::P
    lidx::L
end

function vtu_output(f, u::AbstractVecOrMat, ustr::AbstractString, cache::VTUStructuredCellDataCache)
    vtkfile = vtk_grid(f, cache.pts, cache.cells)
    map(i -> cache.dat[i] = u[cache.tagmap[cache.etag[i]]], cache.lidx)
    vtk_cell_data(vtkfile, cache.dat, ustr)
    vtk_save(vtkfile)
    return vtkfile
end

function vtu_output(f, u::AbstractVector{<:AbstractVecOrMat}, ustr::AbstractVector{<:AbstractString}, cache::VTUStructuredCellDataCache)
    vtkfile = vtk_grid(f, cache.pts, cache.cells)
    for (_u, _ustr) in zip(u, ustr)
        map(i -> cache.dat[i] = _u[cache.tagmap[cache.etag[i]]], cache.lidx)
        vtk_cell_data(vtkfile, cache.dat, _ustr)
    end
    vtk_save(vtkfile)
    return vtkfile
end

function vtu_output(f, t::AbstractVector, u::AbstractArray{<:Number}, ustr::AbstractString, cache::VTUStructuredCellDataCache)
    fmt = "%0$(ndigits(length(t)))d"
    paraview_collection(f) do pvd
        for i = 1: length(t)
            uslice = selectdim(u, ndims(u), i)
            vtkfile = vtu_output(f * sprintf1(fmt, i), uslice, ustr, cache)
            collection_add_timestep(pvd, vtkfile, t[i])
        end
    end
end

function vtu_output(f, t::AbstractVector, u::AbstractVector{<:AbstractArray}, ustr::AbstractVector{<:AbstractString}, cache::VTUStructuredCellDataCache)
    fmt = "%0$(ndigits(length(t)))d"
    paraview_collection(f) do pvd
        for i = 1: length(t)
            us = [selectdim(_u, ndims(_u), i) for _u in u]
            vtkfile = vtu_output(f * sprintf1(fmt, i), us, ustr, cache)
            collection_add_timestep(pvd, vtkfile, t[i])
        end
    end
end
