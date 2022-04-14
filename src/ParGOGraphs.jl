module ParGOGraphs

using DataStructures: Queue, enqueue!, dequeue!, PriorityQueue, isempty

import Base: show, getindex

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
include("output/dot.jl")
include("utils.jl")

# Extensions of the show() function
function show(io::IO, edge::Edge{T}) where {T<:Integer}
    return print(io, "Edge($(source(edge)), $(destination(edge)), $(weight(edge)))")
end

function show(io::IO, G::Union{Graph{T},WeightedGraph{T}}) where {T<:Integer}
    return print(io, "{$(nv(G)), $(ne(G))} undirected graph")
end

function show(io::IO, G::Union{DiGraph{T},WeightedDiGraph{T}}) where {T<:Integer}
    return print(io, "{$(nv(G)), $(ne(G))} directed graph")
end

# The code below is taken from JuMP.jl to export everything, 
# you can see more about it in https://github.com/jump-dev/JuMP.jl/blob/master/src/JuMP.jl.
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__, all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end

end
