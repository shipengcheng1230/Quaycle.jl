export vtu_output

const gmshcelltype2vtkcelltype = Dict(
    1 => VTKCellTypes.VTK_LINE,
    3 => VTKCellTypes.VTK_QUAD,
    4 => VTKCellTypes.VTK_TETRA,
    5 => VTKCellTypes.VTK_HEXAHEDRON,
)

abstract type ParaviewOutputCache end

struct VTKStructuredScalarConversionCache{D, V, A, C, P, L} <: ParaviewOutputCache
    tagmap::D
    etag::V
    dat::A
    cells::C
    pts::P
    lidx::L
end

struct VTKUnStructuredCache{C, P} <: ParaviewOutputCache
    cells::C
    pts::P
end

_map_data(cache::VTKStructuredScalarConversionCache, u::AbstractArray) = map(i -> cache.dat[i] = u[cache.tagmap[cache.etag[i]]], cache.lidx)

function _write_cell_data(vtkfile, u, ustr, cache::VTKStructuredScalarConversionCache)
    _map_data(cache, u)
    vtk_cell_data(vtkfile, cache.dat, ustr)
end

_write_cell_data(vtkfile, u, ustr, cache::VTKUnStructuredCache) = vtk_cell_data(vtkfile, u, ustr)

function vtk_output(f, u::AbstractVector{<:AbstractVecOrMat}, ustr::AbstractVector{<:AbstractString}, cache::ParaviewOutputCache)
    vtk_grid(f, cache.pts, cache.cells) do vtk
        for (_u, _ustr) in zip(u, ustr)
            _write_cell_data(vtk, _u, _ustr, cache)
        end
    end
end

function vtk_output(f, t::AbstractVector, u::AbstractVector{<:AbstractArray}, ustr::AbstractVector{<:AbstractString}, cache::ParaviewOutputCache)
    fmt = "%0$(ndigits(length(t)))d"
    paraview_collection(f) do pvd
        for i = 1: length(t)
            us = [selectdim(_u, ndims(_u), i) for _u in u]
            vtk_grid(f * sprintf1(fmt, i), cache.pts, cache.cells) do vtk
                for (_u, _ustr) in zip(us, ustr)
                    _write_cell_data(vtk, _u, _ustr, cache)
                end
                collection_add_timestep(pvd, vtk, t[i])
            end
        end
    end
end
