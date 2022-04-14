"""
    _update_weights(G::Union{AbstractWeightedGraph, AbstractWeightedDiGraph}, state::MaxFlowState)

Update the weights of the graph based on the flow in given in `state`.

# Arguments
- `G`: Graph to be updated.
- `state`: `MaxFlowState` after some iterations of a maximum flow algorithm.
"""
function _update_weights(
    G::Union{AbstractWeightedGraph,AbstractWeightedDiGraph},
    state::MaxFlowState,
)
    @inbounds for i in eachindex(G.adjacency_list)
        G.weights[i] .-= state.flows[i]
    end
end

"""
    minimum_cut!(G::Union{AbstractWeightedGraph, AbstractWeightedDiGraph}, s::Int, t::Int)::Tuple

Calculates the minimum cut in a graph `G` computing its maximum flow, using Edmonds-Karp algorithm [1] based on the Max-Flow Min-Cut Theorem [2]. The algorithm uses BFS to look for paths. This function changes the weights of the graph `G` such that it corresponds to the updated capacities after the final iteration of the Edmonds-Karp algorithm.

[1] https://en.wikipedia.org/wiki/Edmonds%E2%80%93Karp_algorithm.
[2] https://en.wikipedia.org/wiki/Max-flow_min-cut_theorem.
"""
function minimum_cut!(
    G::Union{AbstractWeightedGraph,AbstractWeightedDiGraph},
    s::Int,
    t::Int,
)::Tuple
    state = ford_fulkerson(G, s, t)

    _update_weights(G, state)

    S_1 = breadth_first_search(G, s, true)
    S_2 = Set(setdiff(1:nv(G), S_1))

    return S_1, S_2
end

"""
    minimum_cut(G::Union{AbstractWeightedGraph, AbstractWeightedDiGraph}, s::Int, t::Int)::Tuple

Calculates the minimum cut in a graph `G` computing its maximum flow, using Edmonds-Karp algorithm [1] based on the Max-Flow Min-Cut Theorem [2]. The algorithm uses BFS to look for paths.

[1] https://en.wikipedia.org/wiki/Edmonds%E2%80%93Karp_algorithm.
[2] https://en.wikipedia.org/wiki/Max-flow_min-cut_theorem.

# Arguments
- `G`: Graph from where the minimum cut is going to be obtained.
- `s`: The source vertex.
- `t`: The sink vertex.
"""
function minimum_cut(G::Union{AbstractWeightedGraph,AbstractWeightedDiGraph}, s::Int, t::Int)::Tuple
    state = ford_fulkerson(G, s, t)
    copy_G = deepcopy(G)

    _update_weights(copy_G, state)

    S_1 = breadth_first_search(copy_G, s, true)
    S_2 = Set(setdiff(1:nv(copy_G), S_1))

    return S_1, S_2
end

"""
    minimum_cut(G::Union{AbstractWeightedGraph, AbstractWeightedDiGraph}, state::MaxFlowState, s::Int)::Tuple

Calculates the minimum cut from a given MaxFlowState.

# Arguments
- `state`: A predefined maximum flow state.
"""
function minimum_cut(
    G::Union{AbstractWeightedGraph,AbstractWeightedDiGraph},
    state::MaxFlowState,
    s::Int,
)::Tuple
    copy_G = deepcopy(G)

    _update_weights(copy_G, state)

    S_1 = breadth_first_search(copy_G, s, true)
    S_2 = Set(setdiff(1:nv(copy_G), S_1))

    return S_1, S_2
end
