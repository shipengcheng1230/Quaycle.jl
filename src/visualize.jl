## These functions are not fully tested, use with caution.
export vtk_output

const gmshcelltype2vtkcelltype = Dict(
    1 => VTKCellTypes.VTK_LINE,
    2 => VTKCellTypes.VTK_TRIANGLE,
    3 => VTKCellTypes.VTK_QUAD,
    4 => VTKCellTypes.VTK_TETRA,
    5 => VTKCellTypes.VTK_HEXAHEDRON,
)

abstract type ParaviewOutputCache end

struct VTKStructuredScalarConversionCache{D, V, A, C, P, L} <: ParaviewOutputCache
    tagmap::D
    tag::V
    dat::A
    cells::C
    pts::P
    lidx::L
end

struct VTKUnStructuredCache{C, P, V} <: ParaviewOutputCache
    cells::C
    pts::P
    tag::V
end

# This works only for scalar data, meaning `u` will be flattened
_map_data(cache::VTKStructuredScalarConversionCache, u::AbstractArray) = map(i -> cache.dat[i] = u[cache.tagmap[cache.tag[i]]], cache.lidx)
_map_data(cache::VTKStructuredScalarConversionCache, u::Number) = cache.dat .= u

function _write_cell_data(vtkfile, u, ustr, cache::VTKStructuredScalarConversionCache)
    _map_data(cache, u)
    vtk_cell_data(vtkfile, cache.dat, ustr)
end

_write_cell_data(vtkfile, u, ustr, cache::VTKUnStructuredCache) = vtk_cell_data(vtkfile, u, ustr)

"""
    vtk_output(f, u::AbstractVector{<:AbstractVecOrMat},
        ustr::AbstractVector{<:AbstractString}, cache::ParaviewOutputCache)

Write results to single-block VTU file.

## Arguments
- `f`: output file name
- `u::AbstractVector`, list of results to be written
- `ustr::AbstractVector{<:AbstractString}`: list results names to be assigned
- `cache`: cache of data conversion, cell information and nodes information
"""
function vtk_output(f, u::AbstractVector, ustr::AbstractVector{<:AbstractString}, cache::ParaviewOutputCache)
    vtk_grid(f, cache.pts, cache.cells) do vtk
        for (_u, _ustr) in zip(u, ustr)
            _write_cell_data(vtk, _u, _ustr, cache)
        end
    end
end

"""
    vtk_output(f, t::AbstractVector, u::AbstractVector{<:AbstractVecOrMat},
        ustr::AbstractVector{<:AbstractString}, cache::ParaviewOutputCache)

Write time-series results to single-block paraview collection file.

## Arguments
- `f`: output file name
- `t::AbstractVector`: time stamp vector
- `u::AbstractVector{<:AbstractArray}`, list of results to be written
- `ustr::AbstractVector{<:AbstractString}`: list results names to be assigned
- `cache`: cache of data conversion, cell information and nodes information
"""
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

"""
    vtk_output(f, u, ustr, cache::AbstractVector{<:ParaviewOutputCache})

Write results to multiple-block VTM file.

## Arguments
- `f`: output file name
- `u`, list of block to be written, each contains a list of results
    in that block
- `ustr`: list of block to be assigned, each contains a list of results names
    to data in that block
- `cache`: list cache corresponding to each block
"""
function vtk_output(f, u, ustr, cache::AbstractVector{<:ParaviewOutputCache})
    vtk_multiblock(f) do vtm
        for (_u, _ustr, _cache) in zip(u, ustr, cache)
            vtkfile = vtk_grid(vtm, _cache.pts, _cache.cells)
            for (u′, ustr′) in zip(_u, _ustr)
                _write_cell_data(vtkfile, u′, ustr′, _cache)
            end
        end
    end
end

"""
    vtk_output(f, t, u, ustr, cache::AbstractVector{<:ParaviewOutputCache})

Write results to multiple-block paraview collection file.

## Arguments
- `f`: output file name
- `t::AbstractVector`: time stamp vector
- `u`, list of different block data to be written, each contains a list of results
    in that block
- `ustr`: list of block names to be assigned, each contains a list of results names
    to data in that block
- `cache`: list cache corresponding to each block
"""
function vtk_output(f, t, u, ustr, cache::AbstractVector{<:ParaviewOutputCache})
    fmt = "%0$(ndigits(length(t)))d"
    paraview_collection(f) do pvd
        for i = 1: length(t)
            vtk_multiblock(f * sprintf1(fmt, i)) do vtm
                for (_u, _ustr, _cache) in zip(u, ustr, cache)
                    u′ = [selectdim(x, ndims(x), i) for x in _u]
                    vtkfile = vtk_grid(vtm, _cache.pts, _cache.cells)
                    for (u′′, ustr′′) in zip(u′, _ustr)
                        _write_cell_data(vtkfile, u′′, ustr′′, _cache)
                    end
                end
                collection_add_timestep(pvd, vtm, t[i])
            end
        end
    end
end

"""
    vtk_output(f, p::AbstractProperty, cache::ParaviewOutputCache)

Store simulation property, including all its fields, to VTU file.

## Arguments
- `f`: output file name
- `p`: property struct
- `cache`: cache of data conversion, cell information and nodes information

## Notice
Not all `AbstractProperty` are supported. Only those with well defined `fieldnames`
    are. For multi-block domain property, save them separately.
"""
vtk_output(f, p::AbstractProperty, cache::ParaviewOutputCache) = vtk_output(f, _get_field_and_name(p)..., cache)

function _get_field_and_name(p::AbstractProperty)
    allname = collect(fieldnames(p))
    allu = [getfield(p, Symbol(pn)) for pn in allname]
    return allu, allname
end
