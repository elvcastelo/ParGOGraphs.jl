"""
    DijkstraState

Structure representing the state of the Dijkstra search algorithm.
"""
mutable struct DijkstraState
    predecessors::Vector{Int}
    distances::Vector{Float64}
    path_length::Float64
end

"""
    DijkstraState(G::AbstractGraph, s::Int)

Initialize the Dijkstra state in the graph `G` starting from node `s`.
"""
function DijkstraState(G::AbstractGraph, s::Int)
    N = nv(G)
    predecessors = zeros(Int, N)
    distances = [Inf for _ in 1:N]
    distances[s] = 0
    state = DijkstraState(predecessors, distances, 0.0)

    return state
end

"""
    dijkstra(G::Union{AbstractWeightedDiGraph, AbstractWeightedGraph}, s::Int, t::Int)::DijkstraState

Implements Dijkstra's algorithm [1] to find the shortest (`s`, `t`)-path in a given weighted graph `G`. The function returns a `DijkstraState` containing information about the distances to the nodes, the shortest path cost and information about every vertex predecessor.

# References

[1] T. H. Cormen, Ed., Introduction to algorithms, 3rd ed. Cambridge, Mass: MIT Press, 2009.

See also: [`get_path`](@ref).
"""
function dijkstra(
    G::Union{AbstractWeightedDiGraph,AbstractWeightedGraph},
    s::Int,
    t::Int;
    ignore_zeros::Bool = false,
)::DijkstraState
    state = DijkstraState(G, s)
    # S = Set{Int}()
    Q = PriorityQueue{Int,Float64}()

    # Populate Q such that Q = G.V
    @inbounds for i in 1:nv(G)
        enqueue!(Q, i, state.distances[i])
    end

    while !isempty(Q)
        u = dequeue!(Q)
        # push!(S, u)
        adjlist = G.adjacency_list[u]
        @inbounds for i in eachindex(adjlist)
            v = adjlist[i]
            weight = G.weights[u][i]

            (weight == 0 && ignore_zeros) && continue

            if state.distances[v] > state.distances[u] + weight
                state.distances[v] = state.distances[u] + weight
                state.predecessors[v] = u
                Q[v] = state.distances[v]
            end
        end
    end

    state.path_length = state.distances[t]

    return state
end

"""
    get_path(state::DijkstraState, dst::Int)::Vector{Tuple}

Obtain the shortest (s, `dst`)-path from the `state` returned after running the dijkstra algorithm.

See also: [`dijkstra`](@ref).
"""
function get_path(state::DijkstraState, dst::Int)::Vector{Tuple{Int,Int}}
    state_path = Vector{Tuple}()
    v = dst
    u = state.predecessors[dst]

    while u != 0
        push!(state_path, (u, v))
        v = u
        u = state.predecessors[u]
    end

    return state_path
end
