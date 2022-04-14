"""
    graphplot(G::AbstractGraph; kwargs...)

Expands `graphplot` to work with PargoGraphs.AbstractGraph, using its adjacency matrix representation.

The recognized keyword arguments are the same from `graphplot`.

See also: [`adjacency_matrix`](@ref).
"""
function graphplot(G::AbstractGraph; kwargs...)
    return graphplot(adjacency_matrix(G); kwargs...)
end
