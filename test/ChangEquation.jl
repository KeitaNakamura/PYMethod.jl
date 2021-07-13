@testset "Chang's equation" begin
    E = 2e8
    I = 0.0002507
    calc_moment(eq, z) = -ForwardDiff.derivative(z′ -> ForwardDiff.derivative(z′′ -> PYMethod.calculate_deflection(eq, z′′), z′), z) * (E*I)
    calc_shearforce(eq, z) = ForwardDiff.derivative(z′ -> calc_moment(eq, z′), z) * (-1) # the direction of z axis is downward
    for eq in (ChangEquation(0, 20; z_0 = 19, F_t = 10, D = 0.6, E, I, k = 3750),
               ChangEquation(0, 20; z_0 = 19, F_t = 10, D = 0.6, E, I, k = 3750, head_free = false),)
        if eq.M_t == 0
            @test eq.θ_t != 0
        else
            @test eq.θ_t == 0
        end
        for z in 0:0.1:20
            @test PYMethod.calculate_moment(eq, z) ≈ calc_moment(eq, z)
            @test PYMethod.calculate_shearforce(eq, z) ≈ calc_shearforce(eq, z)
        end
        # check continuity
        z_0₊ = eq.z_0 + eps(eq.z_0)
        z_0₋ = eq.z_0 - eps(eq.z_0)
        @test PYMethod.calculate_deflection(eq, z_0₊) ≈ PYMethod.calculate_deflection(eq, z_0₋)
        @test PYMethod.calculate_moment(eq, z_0₊) ≈ PYMethod.calculate_moment(eq, z_0₋)
        @test PYMethod.calculate_shearforce(eq, z_0₊) ≈ PYMethod.calculate_shearforce(eq, z_0₋)
    end
end
