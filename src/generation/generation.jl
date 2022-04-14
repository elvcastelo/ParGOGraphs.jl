"""
    erdos_renyi(nvertices::Int, nedges::Int; kwargs...)::AbstractGraph

Creates a random graph with `nvertices` vertices and `nedges` edges using the Erdős–Rényi model [1]. By default the generated graph is undirected. The generation follows an uniform distribution, both for choice of edges and weights.

The recognized keyword arguments are the following:

- `is_directed`: Indicates the graph being generated is directed.
- `is_weighted`: Indicates the graph being generated is weighted.
- `limit_degree`: Gives a upper limit for the degrees of each vertex in the graph.
- `lbound`: Gives a lower bound for the graph weights.
- `ubound`: Gives a upper bound for the graph weights.

# References

[1] https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model

See also [`layered_graph`, `disk_intersection_graph`](@ref).
"""
function erdos_renyi(
    nvertices::Int,
    nedges::Int;
    is_directed::Bool = false,
    is_weighted::Bool = false,
    limit_degree::Int = nvertices,
    lbound::Float64 = 1.0,
    ubound::Float64 = 100.0,
)::AbstractGraph
    total_edges = nedges

    if !is_directed
        total_edges += nedges
    end

    num = nvertices * (nvertices - 1)
    max_edges = is_directed ? num : num / 2

    if total_edges > max_edges
        max_edges_string = is_directed ? "|V| × (|V| - 1)" : "(|V| × (|V| - 1)) / 2"
        graph_type = is_directed ? "digraph" : "graph"

        throw(ErrorException("there can't be more than $max_edges_string edges in a $graph_type"))
    end

    if is_weighted
        if lbound > ubound
            throw(ErrorException("the lower bound is greater than the upper bound"))
        end

        return _weighted_erdos_renyi(
            nvertices,
            nedges,
            lbound,
            ubound,
            is_directed,
            limit_degree,
        )
    end

    return _erdos_renyi(nvertices, nedges, is_directed)
end

"""
    _erdos_renyi(nvertices::Int, nedges::Int[, is_directed::Bool=false])::AbstractGraph

Create a random non-weighted simple (di)graph with `nvertices` vertices and `nedges` edges using the Erdős–Rényi model. Optionally, a digraph can be created by specifying `is_direct = true`.
"""
function _erdos_renyi(nvertices::Int, nedges::Int, is_directed::Bool = false)::AbstractGraph
    G = is_directed ? DiGraph(nvertices) : Graph(nvertices)

    # If G is an undirected graph then we are going to add two edges (i, j) and (j, i)
    max_edges = is_directed ? nedges : 2 * nedges

    while ne(G) < max_edges
        source = rand(1:nvertices)
        destination = rand(1:nvertices)

        if source != destination && destination ∉ G.adjacency_list[source]
            add_edge!(G, source, destination)
        end
    end

    return G
end

"""
    _weighted_erdos_renyi(nvertices::Int, nedges::Int[, lbound::Float64=1., ubound::Float64=100., is_directed::Bool=false]; kwargs...)::AbstractGraph

Create a random weighted simple (di)graph with `nvertices` vertices and `nedges` edges using the Erdős–Rényi model. Optionally, a digraph can be created by specifying `is_direct = true`. The lower bound of the edges weights is given by `lbound` with default value of `1.0`, the upper bound is defined by `ubound` with default value `100.0`.

The recognized keyword arguments are the following:

- `is_direct`: Creates a digraph.
- `limit_degree`: Gives a upper limit for the degrees of each vertex in the graph.
"""
function _weighted_erdos_renyi(
    nvertices::Int,
    nedges::Int,
    lbound::Float64,
    ubound::Float64,
    is_directed::Bool,
    limit_degree::Int,
)::AbstractGraph
    G = is_directed ? WeightedDiGraph(nvertices) : WeightedGraph(nvertices)

    max_edges = is_directed ? nedges : 2 * nedges

    while ne(G) < max_edges
        source = rand(1:nvertices)
        destination = rand(1:nvertices)

        if source != destination &&
           destination ∉ G.adjacency_list[source] &&
           length(G.adjacency_list[source]) < limit_degree
            weight = rand(lbound:ubound)
            add_edge!(G, source, destination, weight)
        end
    end

    return G
end

"""
    _weighted_layered_graph(n::Int, m::Int, is_directed::Bool, lbound::Float64, ubound::Float64):AbstractGraph

Generates a weighted layered graph with `n` columns and `m` vertices per column such that each weight `w` is in the interval [`lbound`, `ubound`].
"""
function _weighted_layered_graph(
    n::Int,
    m::Int,
    is_directed::Bool,
    lbound::Float64,
    ubound::Float64,
)
    :AbstractGraph
    nv = 2 + (n * m)
    G = is_directed ? WeightedDiGraph(nv) : WeightedGraph(nv)

    for j in 2:m+1
        weight = rand(lbound:ubound)
        add_edge!(G, 1, j, weight)
    end

    for j in 1:m
        weight = rand(lbound:ubound)
        v = j + (m * (n - 1)) + 1
        add_edge!(G, v, nv, weight)
    end

    @inbounds for i in 1:n-1
        for j in 1:m
            v = j + (m * (i - 1)) + 1
            for k in 1:m
                weight = rand(lbound:ubound)
                u = k + (m * i) + 1
                add_edge!(G, v, u, weight)
            end
        end
    end

    return G
end

"""
    _layered_graph(n::Int, m::Int; is_directed::Bool)::AbstractGraph

Generates a non-weighted layered graph with `n` columns and `m` vertices per column.
"""
function _layered_graph(n::Int, m::Int, is_directed::Bool)::AbstractGraph
    nv = 2 + (n * m)
    G = is_directed ? DiGraph(nv) : Graph(nv)

    @inbounds for j in 2:m+1
        add_edge!(G, 1, j)
    end

    @inbounds for j in 1:m
        v = j + (m * (n - 1)) + 1
        add_edge!(G, v, nv)
    end

    @inbounds for i in 1:n-1
        for j in 1:m
            v = j + (m * (i - 1)) + 1
            for k in 1:m
                u = k + (m * i) + 1
                add_edge!(G, v, u)
            end
        end
    end

    return G
end

"""
    layered_graph(n::Int, m::Int; is_directed::Bool = false, is_weighted::Bool = true)::AbstractGraph

Generates a layered graph with `n` columns and `m` vertices per column. The generated graph can be directed by setting `is_directed = true` and weighted with `is_weighted = true`. 
"""
function layered_graph(
    n::Int,
    m::Int;
    is_directed::Bool = false,
    is_weighted::Bool = true,
    lbound::Float64 = 1.0,
    ubound::Float64 = 100.0,
)::AbstractGraph
    if is_weighted
        return _weighted_layered_graph(
            n,
            m,
            is_directed,
            lbound,
            ubound,
        )
    end

    return _layered_graph(n, m, is_directed)
end

"""
    _disk_intersection_graph(centres::Vector{Tuple{Int,Int}})::AbstractGraph

Generate a non weighted unit disk graph given centres `centres`. 
"""
function _disk_intersection_graph(centres::Vector{Tuple{Int,Int}})::AbstractGraph
    n = length(centres)

    G = Graph(n)

    @inbounds for i in eachindex(centres)
        for j in i+1:n
            if norm(centres[i] .- centres[j], Inf) <= 1
                add_edge!(G, i, j)
            end
        end
    end

    return G
end

"""
    _weighted_disk_intersection_graph(centres::Vector{Tuple{Int,Int}}, lbound::Float64, ubound::Float64)

Generate a weighted unit disk graph given centres `centres` with weight lower bound `lbound` and upper bound `ubound`.
"""
function _weighted_disk_intersection_graph(
    centres::Vector{Tuple{Int,Int}},
    lbound::Float64,
    ubound::Float64,
)
    n = length(centres)

    G = WeightedGraph(n)

    @inbounds for i in eachindex(centres)
        for j in i+1:n
            if norm(centres[i] .- centres[j], Inf) <= 1
                weight = rand(lbound:ubound)
                add_edge!(G, i, j, weight)
            end
        end
    end

    return G
end

# TODO: Improve the algorithm complexity by using the work of Bentley et al (1977)
"""
    function disk_intersection_graph(n::Int, b::Int, h::Int; is_weighted::Bool = false, lbound::Float64 = 1.0, ubound::Float64 = 100.0)::AbstractGraph

Generate a unit disk graph [1] with `n` vertices based on a rectangle `b` x `h`. The graph can be weighted by setting `is_weighted = true`. A lower and upper bound can be specified with `lbound` and `ubound` respectively. 

# References
[1] https://en.wikipedia.org/wiki/Unit_disk_graph
"""
function disk_intersection_graph(
    n::Int,
    b::Int,
    h::Int,
    is_weighted::Bool = false,
    lbound::Float64 = 1.0,
    ubound::Float64 = 100.0,
)::AbstractGraph
    centres = [(rand(1:b), rand(1:h)) for _ in 1:n]
    is_weighted && return _weighted_disk_intersection_graph(centres, lbound, ubound)
    return _disk_intersection_graph(centres)
end
