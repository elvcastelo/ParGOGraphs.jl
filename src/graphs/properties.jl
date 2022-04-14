"""
    outneighbors(G::AbstractGraph, vertex::Int)::Vector{Int}

Return the neighbors of the vertex `vertex` in the graph `G`, such that an edge (`vertex`, `neighbor`) for all `neighbor` in `neighbors` exists. This function corresponds to returning the list `G.adjacency_list[vertex]` stored in the graph structure.

See also [`inneighbors`](@ref).
"""
outneighbors(G::AbstractGraph, vertex::Int)::Vector{Int} = return G.adjacency_list[vertex]

"""
    inneighbors(G::AbstractUndirectedGraph, vertex::Int)::Vector{Int}

Return the neighbors of the vertex `vertex` such that an edge (`neighbor`, `vertex`) for all `neighbor` in `neighbors` exists. For an undirected graphs this function behaves the same as `outneighbors`.

See also [`outneighbors`](@ref).
"""
inneighbors(G::AbstractUndirectedGraph, vertex::Int)::Vector{Int} = return G.adjacency_list[vertex]

"""
    inneighbors(G::AbstractDiGraph, vertex::Int)::Vector{Int}

For digraphs this function return `G.back_adjacency_list[vertex]`.
    
See also [`outneighbors`](@ref).
"""
inneighbors(G::AbstractDiGraph, vertex::Int)::Vector{Int} = return G.back_adjacency_list[vertex]

"""
    nv(G::AbstractGraph)::Int

Return the number of vertices of a graph `G`.
"""
nv(G::AbstractGraph)::Int = return G.nvertices

"""
    ne(G::AbstractGraph)::Int

Return the number of edges of a graph `G`.
"""
ne(G::AbstractGraph)::Int = return G.nedges
