module ParGOGraphs

using DataStructures: Queue, enqueue!, dequeue!, PriorityQueue, isempty

import Base: getindex

# TODO: Improve on this
export
    # Base
    AbstractGraph,
    AbstractUndirectedGraph,
    AbstractDiGraph,
    AbstractWeightedGraph,
    AbstractWeightedDiGraph,
    Graph,
    DiGraph,
    WeightedGraph,
    WeightedDiGraph,
    AbstractEdge,
    Edge,
    add_edge!,
    edges,
    has_edge,
    source,
    destination,
    weight,
    getindex,
    outneighbors,
    inneighbors,
    nv,
    ne,

    # Generation
    erdos_renyi,
    layered_graph,
    disk_intersection_graph,

    # Maximum-Flow
    MaxFlowState,
    ford_fulkerson,
    minimum_cut,
    minimum_cut!,

    # Plot
    graphplot,

    # Traversal
    DijkstraState,
    dijkstra,
    breadth_first_search,
    get_path,

    # Utils
    adjacency_matrix

include("graphs/graph.jl")
include("graphs/digraph.jl")
include("graphs/edge.jl")
include("graphs/properties.jl")

include("generation/generation.jl")

include("maximum_flow/ford_fulkerson.jl")
include("maximum_flow/minimum_cut.jl")

include("plot/plot.jl")

include("traversal/breadth.jl")
include("traversal/dijkstra.jl")

include("utils.jl")

end
