"""
    save_dot(G::AbstractGraph, path::String)

Saves the graph `G` to `path` as a dotfile [1].

# References

[1] https://en.wikipedia.org/wiki/DOT_(graph_description_language)
"""
function save_dot(G::AbstractGraph, path::String)
    name = typeof(G) <: AbstractUndirectedGraph ? "graph" : "digraph"
    arrow = typeof(G) <: AbstractUndirectedGraph ? " -- " : " -> "
    edges_list = ""

    for i in 1:nv(G)
        adjlist = G.adjacency_list[i]
        for j in 1:length(adjlist)
            (typeof(G) <: AbstractUndirectedGraph && i > adjlist[j]) && continue

            if typeof(G) <: AbstractWeightedGraph || typeof(G) <: AbstractWeightedDiGraph
                edges_list *=
                    string(i) *
                    arrow *
                    string(adjlist[j]) *
                    " [label=$(round(Int, G.weights[i][j]))] \n\t"
            else
                edges_list *= string(i) * arrow * string(adjlist[j]) * "\n\t"
            end
        end
    end

    dot = """
    $name {
        node[shape=circle]
        $edges_list
    }
    """

    open(path, "w") do io
        return write(io, dot)
    end
end