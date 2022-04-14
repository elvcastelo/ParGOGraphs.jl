"""
    MaxFlowState

A structure to represent the current state of a maximum flow algorithm.
"""
mutable struct MaxFlowState
    flows::Vector{Vector{Float64}}
    capacities::Vector{Vector{Float64}}
    max_flow::Float64
end

"""
    flow(state::MaxFlowState, G::AbstractGraph, s, t)::Float64

Returns the current flow in (`s`, `t`) in the graph `G` from the MaxFlowState `state`.
"""
function flow(state::MaxFlowState, G::AbstractGraph, s, t)::Float64
    index = findall(G.adjacency_list[s] .== t)[1]
    return state.flows[s][index]
end

"""
    capacity(state::MaxFlowState, G::AbstractGraph, s, t)::Float64

Returns the current capacity in the edge `(s, t)`.

# Arguments
- `state`: Current maximum flow state in the graph `G`.
- `G`: Reference graph.
- `s`: Source vertex.
- `t`: Destination vertex.

# Returns
- `capacity`: Current capacity in the edge `(s, t)`.
"""
function capacity(state::MaxFlowState, G::AbstractGraph, s, t)::Float64
    index = findall(G.adjacency_list[s] .== t)
    return state.capacities[s][index...]
end

"""
    flow_capacity(state::MaxFlowState, G::AbstractGraph, s, t)::Tuple{Float64, Float64}

Returns the current flow and capacity in the edge `(s, t)`.

# Arguments
- `state`: Current maximum flow state in the graph `G`.
- `G`: Reference graph.
- `s`: Source vertex.
- `t`: Destination vertex.

# Returns
- `flow`: Current flow in the edge `(s, t)`.
- `capacity`: Current capacity in the edge `(s, t)`.
"""
function flow_capacity(state::MaxFlowState, G::AbstractGraph, s, t)::Tuple{Float64,Float64}
    index = findall(G.adjacency_list[s] .== t)
    return state.flows[s][index][1], state.capacities[s][index][1]
end

"""
    residual_graph(G::AbstractGraph, state::MaxFlowState)::AbstractGraph

Calculate the residual graph of the graph `G`. The residual graph is a weighted graph where its weights represents its capacities, being a more compact representation. This functions assumes the graph has no antiparallel edges [1].

References:

[1] Cormen TH, editor. Introduction to algorithms. 3rd ed. Cambridge, Mass: MIT Press; 2009. 1292 p. 

# Arguments
- `G`: Graph which is to create the residual graph from.
- `state`: Initial maximum flow state.

# Returns
- `G_f`: Residual graph created from `G`.
"""
function residual_graph(G::AbstractDiGraph, state::MaxFlowState)::AbstractGraph
    G_f = WeightedDiGraph(nv(G))
    for edge in edges(G)
        s = source(edge)
        t = destination(edge)

        capacity_st = capacity(state, G, s, t)

        add_edge!(G_f, s, t, capacity_st)
        add_edge!(G_f, t, s, 0.0) # Initially, there is no flow
    end

    return G_f
end

"""
    update_flow!(state::MaxFlowState, G::AbstractGraph, G_f::AbstractGraph, s, t, flow)

Updates the flow in the MaxFlowState and `G_f`, the residual graph using the maximum flow from the augmenting path found in `G_f`. It uses the reference graph `G` to calculate the position for the MaxFlowState.

# Arguments
- `state`: Current maximum flow state.
- `G`: Reference graph.
- `G_f`: Residual graph.
- `s`: Source vertex.
- `t`: Destination vertex.
- `flow`: Maximum flow passing though the augmenting path.
"""
function update_flow!(state::MaxFlowState, G::AbstractGraph, G_f::AbstractGraph, s, t, flow)
    indexf_st = findall(G_f.adjacency_list[s] .== t)
    index_st = findall(G.adjacency_list[s] .== t)

    indexf_ts = findall(G_f.adjacency_list[t] .== s)

    G_f.weights[s][indexf_st] .-= flow
    state.flows[s][index_st] .+= flow

    return G_f.weights[t][indexf_ts] .+= flow
end

"""
    update_maximum_flow!(state::MaxFlowState, G::AbstractGraph, t::Int)

Update the information about the maximum flow of the current Max Flow state. This functions searches for the destination in the adjacency list and sums the flows that reaches it.

# Arguments
- `state`: MaxFlowState generated after the final iteration of a maximum flow algorithm.
- `G`: Graph where the maximum flow were executed.
- `t`: Sink node.
"""
function update_maximum_flow!(state::MaxFlowState, G::AbstractGraph, t::Int)
    max_flow = 0

    @inbounds for i in 1:nv(G)
        index = findall(G.adjacency_list[i] .== t)
        if length(index) != 0
            max_flow += state.flows[i][index][1]
        end
    end

    return state.max_flow = max_flow
end

# TODO: Make it work on graphs with antiparallel edges
"""
    ford_fulkerson(G::Union{AbstractWeightedGraph, AbstractWeightedDiGraph}, s, t)::MaxFlowState

Implements the Ford-Fulkerson algorithm [1, 3], more specifically the Edmonds-Karp algorithm [2] from the initial vertex `s` to the destination vertex `t`. This function assumes the graph has no antiparallel edges [3].

References:

[1] https://en.wikipedia.org/wiki/Ford%E2%80%93Fulkerson_algorithm
[2] https://en.wikipedia.org/wiki/Edmonds%E2%80%93Karp_algorithm
[3] Cormen TH, editor. Introduction to algorithms. 3rd ed. Cambridge, Mass: MIT Press; 2009. 1292 p. 

# Arguments
- `G`: Weighted graph where the weights represents the capacities of the edges.
- `s`: Source vertex.
- `t`: Destination vertex.

# Returns
- `state`: Maximum flow final state.
"""
function ford_fulkerson(G::AbstractWeightedDiGraph, s::Int, t::Int)::MaxFlowState

    # TODO: Work with antiparallel edges
    @inbounds for i in 1:nv(G), j in outneighbors(G, i)
        j_index = findall(G.adjacency_list[i] .== j)[1]

        weight = G.weights[i][j_index]

        if i in outneighbors(G, j) && weight != 0
            println((i, j))
            throw(ErrorException("the graph can't have antiparallel edges"))
        end
    end

    state = MaxFlowState([zeros(Int, length(G.adjacency_list[i])) for i in 1:nv(G)], G.weights, 0)

    G_f = residual_graph(G, state)

    while true
        path = breadth_first_search(G_f, s, t)

        length(path) == 0 && break

        max_flow = minimum(weight.(path))

        for i in eachindex(path)
            edge = path[i]
            update_flow!(state, G, G_f, source(edge), destination(edge), max_flow)
        end
    end

    update_maximum_flow!(state, G, t)

    return state
end
