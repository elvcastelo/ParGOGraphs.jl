"""
    adjacency_matrix(G::AbstractGraph)::Matrix{Int}

Creates and return an adjacency matrix representing the graph `G` If `G` is big enough the matrix may be too big to keep in memory.
"""
function adjacency_matrix(G::AbstractGraph)::Matrix{Int}
    N = nv(G)
    A = zeros(Int, N, N)

    @inbounds for i in 1:N
        outneigh = outneighbors(G, i)
        A[i, outneigh] .= 1

        if typeof(G) <: AbstractUndirectedGraph
            A[outneigh, i] .= 1
        else
            inneigh = inneighbors(G, i)
            A[inneigh, i] .= 1
        end
    end

    return A
end
