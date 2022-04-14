abstract type AbstractEdge end

"""
    Edge{T<:Real} <: AbstractEdge

A structure representing a weighted edge.
"""
struct Edge{T<:Real} <: AbstractEdge
    source::Int
    destination::Int
    weight::T
end

"""
    Edge(source::Int, destination::Int)::Edge

Creates the edge (`source`, `destination`) with unitary weight.
"""
Edge(source::Int, destination::Int)::Edge = Edge(source, destination, 0)

"""
    source(edge::Edge)::Int

Return the source node from the edge `edge`.

See also [`destination`](@ref), [`weight`](@ref).
"""
source(edge::Edge)::Int = return edge.source

"""
    destination(edge::Edge)::Int

Return the desination node from the edge `edge`.

See also [`source`](@ref), [`weight`](@ref).
"""
destination(edge::Edge)::Int = return edge.destination

"""
    weight(edge::Edge)::Float64

Return the weight of the edge `edge`.

See also [`source`](@ref), [`destination`](@ref).
"""
weight(edge::Edge)::Float64 = return edge.weight

"""
    getindex(edge::Edge, vals...)

Construct a 1-d array with elements of Edge obtained by the given indexes.
"""
function getindex(edge::Edge, vals...)
    a = Vector{Float64}(undef, length(vals))
    @inbounds for i in 1:length(vals)
        vals[i] > 3 && throw(ErrorException("Edge only have three fields"))
        # TODO: A better way to do this?
        if vals[i] == 1
            a[i] = source(edge)
        elseif vals[i] == 2
            a[i] = destination(edge)
        else
            a[i] = weight(edge)
        end
    end

    return a
end
