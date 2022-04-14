"""
    AbstractGraph

Abstract type representing a generic graph.
"""
abstract type AbstractGraph end

"""
    AbstractUndirectedGraph

Abstract type representing an undirected graph.
"""
abstract type AbstractUndirectedGraph <: AbstractGraph end

"""
    AbstractWeightedGraph

Abstract type representing an undirected weighted graph.
"""
abstract type AbstractWeightedGraph <: AbstractUndirectedGraph end

"""
    Graph{T <: Integer} <: AbstractUndirectedGraph

A structure representing an undirected graph. This structures keep track of the number of edges `nedges`, the number of vertices `nvertices` and mantain an adjacency list `adjacency_list` to represent the graph.
"""
mutable struct Graph{T<:Integer} <: AbstractUndirectedGraph
    nedges::Int
    nvertices::Int
    adjacency_list::Vector{Vector{T}}
end

"""
    Graph{T}(nvertices::Int)::Graph where T <: Integer

Creates an undirected graph of type `T` with `nvertices` vertices.
"""
function Graph{T}(nvertices::Int)::Graph where {T<:Integer}
    return Graph{T}(0, nvertices, [Vector{T}() for _ in 1:nvertices])
end

"""
    Graph(nvertices::T)::Graph where T <: Integer

Creates an undirected graph with `nvertices` vertices of type `T`.
"""
function Graph(nvertices::T)::Graph where {T<:Integer}
    return Graph{T}(nvertices)
end

"""
    add_edge!(G::Graph, source::Int, destination::Int)

Add an edge (`source`, `destination`) and (`destination`, `source`) to the undirected graph `G`.
"""
function add_edge!(G::Graph, source::Int, destination::Int)
    @inbounds list_out = G.adjacency_list[source]
    @inbounds list_in = G.adjacency_list[destination]

    index = searchsortedfirst(list_out, destination)
    insert!(list_out, index, destination)

    index = searchsortedfirst(list_in, source)
    insert!(list_in, index, source)

    return G.nedges += 2
end

"""
    _incidence_graph!(G::Graph, A::Matrix{T}) where T <: Integer

Populates the graph `G` based on the incidence matrix `A` of type `T`.
"""
function _incidence_graph!(G::Graph, A::Matrix{T}) where {T<:Integer}
    n = size(A, 2)

    @inbounds for j in 1:n
        indexes = findall(A[:, j] .!= 0)
        for index in indexes
            add_edge!(G, index, j)
        end
    end
end

"""
    _adjacency_graph!(G::Graph, A::Matrix{T}) where T <: Integer

Populates the graph `G` based on the adjacency matrix `A` of type `T`.
"""
function _adjacency_graph!(G::Graph, A::Matrix{T}) where {T<:Integer}
    m = size(A, 1)

    @inbounds for j in 1:m
        # Since the graph is undirected, the matrix is symmetric
        for i in 1:j
            if A[i, j] != 0
                add_edge!(G, i, j)
            end
        end
    end
end

"""
    Graph(A::Matrix{T}[, incidence = false])::Graph where T <: Integer

Creates an undirected graph from an adjacency matrix `A` of type `T`. 

The recognized positional arguments are the following:

- `incidence`: Boolean indicating if the matrix `A` is an incidence matrix.
"""
function Graph(A::Matrix{T}, incidence = false)::Graph where {T<:Integer}
    m, n = size(A)
    adjacency_list = [Vector{T}() for _ in 1:m]

    if incidence
        G = Graph(n, m, adjacency_list)
        _incidence_graph!(G, A)
    else
        G = Graph(0, m, adjacency_list)
        _adjacency_graph!(G, A)
    end

    return G
end

"""
    edges(G::AbstractGraph)::Vector{Edge}

Return the edges of the graph `G` as a vector of `Edge`. Since `G` is a non-weighted graph we consider unitary weight for each edge.

See also [`Edge`](@ref).
"""
function edges(G::AbstractGraph)::Vector{Edge}
    edge_list = Vector{Edge}()

    for i in 1:nv(G)
        list = G.adjacency_list[i]
        for j in eachindex(list)
            push!(edge_list, Edge(i, list[j]))
        end
    end

    return edge_list
end

"""
    WeightedGraph{T <: Integer} <: AbstractWeightedGraph

A structure representing a weighted undirected graph. It keeps track of the number of edges `nedges`, vertices `nvertices` and represents the graph using an adjacency list `adjacency_list`. The indexes of the adjacency list are the same for `weights`.
"""
mutable struct WeightedGraph{T<:Integer} <: AbstractWeightedGraph
    nedges::Int
    nvertices::Int
    adjacency_list::Vector{Vector{T}}
    weights::Vector{Vector{Float64}}
end

"""
    WeightedGraph{T}(nvertices::Int)::AbstractWeightedGraph where T <: Integer

Creates a weighted undirected graph of type `T` with `nvertices` vertices.
"""
function WeightedGraph{T}(nvertices::Int)::AbstractWeightedGraph where {T<:Integer}
    return WeightedGraph{T}(
        0,
        nvertices,
        [Vector{T}() for _ in 1:nvertices],
        [Vector{Float64}() for _ in 1:nvertices],
    )
end

"""
    WeightedGraph(nvertices::T)::AbstractWeightedGraph where T <: Integer

Creates a weighted undirected graph with `nvertices` vertices of type `T`.
"""
function WeightedGraph(nvertices::T)::AbstractWeightedGraph where {T<:Integer}
    return WeightedGraph{T}(nvertices)
end

"""
    add_edge!(G::WeightedGraph, source::T, destination::T[, weight::Float64 = 1.]) where T <: Integer

Add an edge (`source`, `destination`) with weight `weight` to a weighted graph `G`. If no value for `weight` is given the edge is assigned unitary weight.

The recognized keyword arguments are the following:

- `return_index`: Return the index of `destination` in `G.adjacency_list[source]`.
- `return_indexes`: Return the index of `destination` in `G.adjacency_list[source]` and the index of `source` in `G.adjacency_list[destination]`.
"""
function add_edge!(
    G::AbstractWeightedGraph,
    source::T,
    destination::T,
    weight::Float64 = 1.0;
    kwargs...,
) where {T<:Integer}
    @inbounds list_out = G.adjacency_list[source]
    @inbounds list_in = G.adjacency_list[destination]

    index_out = searchsortedfirst(list_out, destination)
    insert!(list_out, index_out, destination)
    insert!(G.weights[source], index_out, weight)

    index_in = searchsortedfirst(list_in, source)
    insert!(list_in, index_in, source)
    insert!(G.weights[destination], index_in, weight)

    G.nedges += 2

    if :return_index in keys(kwargs) && kwargs[:return_index]
        return index_out
    end

    if :return_indexes in keys(kwargs) && kwargs[:return_indexes]
        return index_out, index_in
    end
end

"""
    edges(G::AbstractWeightedGraph)::Vector{Edge}

Returns the edges of a weighted graph `G` as a vector of `Edge`.
"""
function edges(G::AbstractWeightedGraph)::Vector{Edge}
    edge_list = Vector{Edge}()
    @inbounds for i in 1:nv(G)
        list = G.adjacency_list[i]
        for j in eachindex(list)
            push!(edge_list, Edge(i, list[j], G.weights[i][j]))
        end
    end

    return edge_list
end

"""
    has_edge(G::AbstractUndirectedGraph, s::T, t::T) where T <: Integer

Check if the edge (`s`, `t`) exists in the graph `G`.
"""
function has_edge(G::AbstractUndirectedGraph, s::T, t::T)::Bool where {T<:Integer}
    s > nv(G) && throw(ErrorException("the node $s is not part of the graph G"))
    adjlist = G.adjacency_list[s]
    return insorted(t, adjlist)
end
