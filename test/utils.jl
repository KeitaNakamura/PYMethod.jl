@testset "Utilities" begin
    @test PYMethod.second_moment_of_area_closed_ended_pile(1.5) ≈ π*1.5^4/64
    @test PYMethod.second_moment_of_area_open_ended_pile(1.5, 1.2) ≈ π*(1.5^4-1.2^4)/64
    @test_throws Exception PYMethod.second_moment_of_area_open_ended_pile(1.2, 1.5)
end
