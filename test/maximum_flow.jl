using Test
using ParGOGraphs

@testset "Random instance of Maximum Flow" begin
    G = WeightedDiGraph(6)
    add_edge!(G, 1, 2, 16.0)
    add_edge!(G, 1, 3, 13.0)
    add_edge!(G, 2, 4, 12.0)
    add_edge!(G, 3, 2, 4.0)
    add_edge!(G, 3, 5, 14.0)
    add_edge!(G, 4, 3, 9.0)
    add_edge!(G, 4, 6, 20.0)
    add_edge!(G, 5, 4, 7.0)
    add_edge!(G, 5, 6, 4.0)

    @testset "Solve instance" begin
        state = ford_fulkerson(G, 1, 6)

        @test state.max_flow == 23
    end

    @testset "Throw antiparallel edges error" begin
        add_edge!(G, 2, 3, 5.0)

        @test_throws ErrorException ford_fulkerson(G, 1, 6)
    end
end
