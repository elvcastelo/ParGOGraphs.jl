"""
    AbstractDiGraph

An abstract type representing a digraph.
"""
abstract type AbstractDiGraph <: AbstractGraph end

"""
    AbstractWeightedDiGraph

An abstract type representing a weighted digraph.
"""
abstract type AbstractWeightedDiGraph <: AbstractDiGraph end

"""
    DiGraph{T <: Integer} <: AbstractDiGraph

A structure representing a digraph. It keeps track of the number of edges given by `nedges` and vertices, given by `nvertices`. It uses an adjacency list [1] instead of an adjacency matrix to keep track of the conections to optimize memory usage. It also keep its in-neighborhood and out-neighborhood into separate adjacency lists, namely `back_adjacency_list` and `adjacency_list`.

# References

[1] T. H. Cormen, Ed., Introduction to algorithms, 3rd ed. Cambridge, Mass: MIT Press, 2009.
"""
mutable struct DiGraph{T<:Integer} <: AbstractDiGraph
    nedges::Int
    nvertices::Int
    adjacency_list::Vector{Vector{T}}
    back_adjacency_list::Vector{Vector{T}}
end

"""
    DiGraph{T}(nvertices::Int)::AbstractDiGraph where T <: Integer

Creates a digraph of type `T` with `nvertices` vertices.
"""
function DiGraph{T}(nvertices::Int)::AbstractDiGraph where {T<:Integer}
    return DiGraph{T}(
        0,
        nvertices,
        [Vector{T}() for _ in 1:nvertices],
        [Vector{T}() for _ in 1:nvertices],
    )
end

"""
    DiGraph(nvertices::T)::AbstractDiGraph where T <: Integer

Creates a digraph with `nvertices` vertices of type `T`.
"""
function DiGraph(nvertices::T)::AbstractDiGraph where {T<:Integer}
    return DiGraph{T}(nvertices)
end

"""
    WeightedDiGraph{T <: Integer} <: AbstractWeightedDiGraph

A structure representing a weighted DiGraph. It works as the DiGraph struct with an extra field `weights` to store weights of the edges of the graph. This new field follows the same index as `adjacency_list`.
"""
mutable struct WeightedDiGraph{T<:Integer} <: AbstractWeightedDiGraph
    nedges::Int
    nvertices::Int
    adjacency_list::Vector{Vector{T}}
    back_adjacency_list::Vector{Vector{T}}
    weights::Vector{Vector{Float64}}
end

"""
    WeightedDiGraph{T}(nvertices::Int)::AbstractWeightedDiGraph where T <: Integer

Creates a weighted digraph of type `T` with `nvertices` vertices.
"""
function WeightedDiGraph{T}(nvertices::Int)::AbstractWeightedDiGraph where {T<:Integer}
    adjlist = [Vector{T}() for _ in 1:nvertices]
    back_adjlist = [Vector{T}() for _ in 1:nvertices]
    weights = [Vector{Float64}() for _ in 1:nvertices]
    return WeightedDiGraph{T}(0, nvertices, adjlist, back_adjlist, weights)
end

"""
    WeightedDiGraph{T}(nvertices::Int)::AbstractWeightedDiGraph where T <: Integer

Creates a weighted digraph with `nvertices` vertices of type `T`.
"""
function WeightedDiGraph(nvertices::T)::AbstractWeightedDiGraph where {T<:Integer}
    return WeightedDiGraph{T}(nvertices)
end

"""
    add_edge!(G::DiGraph, source::Int, destination::Int)

Add a new edge (`source`, `destination`) to the digraph `G`. The new edge will be inserted into the field `adjacency_list` and `back_adjacency_list`.
"""
function add_edge!(G::AbstractDiGraph, source::Int, destination::Int)
    @inbounds list_out = G.adjacency_list[source]
    @inbounds list_in = G.back_adjacency_list[destination]

    index = searchsortedfirst(list_out, destination)
    insert!(list_out, index, destination)

    index = searchsortedfirst(list_in, source)
    insert!(list_in, index, source)

    return G.nedges += 1
end

"""
    add_edge!(G::AbstractWeightedDiGraph, source::Int, destination::Int, weight::Float64 = 1.; kwargs...)

Add a new edge (`source`, `destination`) with weight `weight` to the weighted digraph `G`. If `weight` is not given we assume it is equal one. 

The recognized keyword arguments are the following:

- `return_index`: Allows to return the index of `destination` in `G.adjacency_list[source]`.
"""
function add_edge!(
    G::AbstractWeightedDiGraph,
    source::Int,
    destination::Int,
    weight::Float64 = 1.0;
    kwargs...,
)
    @inbounds list_out = G.adjacency_list[source]
    @inbounds list_in = G.back_adjacency_list[destination]
    @inbounds list_weights = G.weights[source]

    index_out = searchsortedfirst(list_out, destination)
    insert!(list_out, index_out, destination)
    insert!(list_weights, index_out, weight)

    index_in = searchsortedfirst(list_in, source)
    insert!(list_in, index_in, source)

    G.nedges += 1

    # Return the index of the destination
    if :return_index in keys(kwargs) && kwargs[:return_index]
        return index_out
    end
end

"""
    has_edge(G::AbstractDiGraph, s::T, t::T)::Bool where T <: Integer

Check if the edge (`s`, `t`) exists in the digraph `G`.
"""
function has_edge(G::AbstractDiGraph, s::T, t::T)::Bool where {T<:Integer}
    s > nv(G) && throw(ErrorException("the node $s is not part of the digraph G"))
    adjlist = G.adjacency_list[s]
    return insorted(t, adjlist)
end
