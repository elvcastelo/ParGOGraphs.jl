"""
    _get_path(G::Union{AbstractWeightedGraph, AbstractWeightedDiGraph}, path_state::Vector{Int}, t::Int)::Vector{Edge}

Obtains the path in the weighted graph `G` with destination node `t` given by the `path_state` returned by a traversal algorithm such as Breadth-First Search (BFS) and Depth-First Search (DFS), returning a vector of `Edge`.
"""
function _get_path(
    G::Union{AbstractWeightedGraph,AbstractWeightedDiGraph},
    path_state::Vector{Int},
    t::Int,
)::Vector{Edge}
    p = Vector{Edge}()
    predecessor = path_state[t]

    # If the predecessor is 0 then t is never reached
    predecessor == 0 && return p

    # Obtain the weight of the edge (j, t)
    t_index = findall(G.adjacency_list[predecessor] .== t)[1]
    weight = G.weights[predecessor][t_index]

    push!(p, Edge(predecessor, t, weight))

    @inbounds while predecessor > 0 && path_state[predecessor] != -1
        # Obtain the new predecessor
        j = path_state[predecessor]
        # Obtain the predecessor weight
        weights_list = G.weights[j]
        # Obtain the current predecessor index
        j_index = findall(G.adjacency_list[j] .== predecessor)[1]
        # Obtain the weight of the edge (j, i)
        weight = weights_list[j_index]

        insert!(p, 1, Edge(j, predecessor, weight))

        predecessor = j
    end

    return p
end

"""
    _get_path(path_state::Vector{Int}, t::Int)::Vector{Edge}

Obtains the path from a unitary weighted graph.
"""
function _get_path(path_state::Vector{Int}, t::Int)::Vector{Edge}
    p = Vector{Edge}()
    predecessor = path_state[t]

    # If the predecessor is 0 then t is never reached
    predecessor != 0 && push!(p, Edge(predecessor, t))

    @inbounds while predecessor > 0 && path_state[predecessor] != -1
        j = path_state[predecessor]
        insert!(p, 1, Edge(j, predecessor))
        predecessor = j
    end

    return p
end

"""
    breadth_first_search(G::AbstractGraph, s::Int, t::Int[, ignore::Bool=false])::Vector{Edge}

Implement a Breadth-First Search (BFS) in a graph `G` starting in node `s` with destination `t`. When `G` isn't a weighted graph (corresponding to having edges weights equals to 1) the BFS result corresponds to the shortest s-t-path. Sometimes it is desirable to ignore edges that have weight equal zero (e.g. in maximum flow algorithms when searching in the residual graph), this behavior is controlled by the `ignore` boolean variable.
"""
function breadth_first_search(G::AbstractGraph, s::Int, t::Int, ignore::Bool = true)::Vector{Edge}
    path_state = zeros(Int, nv(G))
    path_state[s] = -1
    visited = Set{Int}()
    Q = Queue{Int}()
    enqueue!(Q, s)

    while !isempty(Q)
        u = dequeue!(Q)
        if u ∉ visited
            neighbors = outneighbors(G, u)

            if ignore
                adjlist = G.adjacency_list[u]
                list_weights = G.weights[u]
            end

            @inbounds for i in eachindex(neighbors)
                if ignore
                    # Get the index of the weight for the neighbor `i`
                    weight_index = adjlist .== neighbors[i]
                    if path_state[neighbors[i]] == 0 && list_weights[weight_index][1] != 0
                        path_state[neighbors[i]] = u
                        enqueue!(Q, neighbors[i])
                    end
                elseif path_state[neighbors[i]] == 0
                    path_state[neighbors[i]] = u
                    enqueue!(Q, neighbors[i])
                end
            end

            push!(visited, u)
        end
    end

    (typeof(G) <: AbstractWeightedGraph || typeof(G) <: AbstractWeightedDiGraph) &&
        return _get_path(G, path_state, t)

    return _get_path(path_state, t)
end

"""
    breadth_first_search(G::AbstractGraph, s::Int, ignore::Bool=false)

Search in the whole graph `G`, starting from node `s`, instead of returning when a certain vertex is found. The algorithm returns the vertices visited instead of a path.
"""
function breadth_first_search(G::AbstractGraph, s::Int, ignore::Bool = false)
    path_state = zeros(Int, nv(G))
    path_state[s] = -1
    visited = Set{Int}()
    Q = Queue{Int}()
    enqueue!(Q, s)

    while !isempty(Q)
        u = dequeue!(Q)
        if u ∉ visited
            neighbors = outneighbors(G, u)

            if ignore
                adjlist = G.adjacency_list[u]
                list_weights = G.weights[u]
            end

            @inbounds for i in eachindex(neighbors)
                if ignore
                    # Get the index of the weight for the neighbor `i`
                    weight_index = adjlist .== neighbors[i]
                    if path_state[neighbors[i]] == 0 && list_weights[weight_index][1] != 0
                        path_state[neighbors[i]] = u
                        enqueue!(Q, neighbors[i])
                    end
                elseif path_state[neighbors[i]] == 0
                    path_state[neighbors[i]] = u
                    enqueue!(Q, neighbors[i])
                end
            end

            push!(visited, u)
        end
    end

    return visited
end
