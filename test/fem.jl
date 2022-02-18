function deflection_chang(H::Real, D::Real, E::Real, I::Real, k_h::Real, h::Real, x::Real)
    β = (k_h*D/(4E*I))^(1/4)
    if x < 0
        H / (6E*I*β^3) * (β^3*x^3 + 3β^3*h*x^2 - 3(1+2β*h)*β*x + 3(1+β*h))
    else
        H / (2E*I*β^3) * exp(-β*x) * ((1+β*h)*cos(β*x) - β*h*sin(β*x))
    end
end

@testset "FEM: Linear spring" begin
    # compare with Chang's equation
    pile = FEPileModel(0, 20, 400)
    D = 0.6
    E = 2e8
    I = 0.0002507
    k_h = 3750
    H = 10
    pile.D .= D
    pile.E .= E
    pile.I .= I
    pile.pycurves .= pycurve(y, z) = z > 19 ? 0 : k_h*y;
    pile.Fext[1] = H
    solve!(pile)
    deflection(z) = deflection_chang(H, D, E, I, k_h, 1, 19-z)
    for i in 1:length(pile.u)
        @test pile.u[i] ≈ deflection(pile.coordinates[i])  atol = 1e-4
    end
end

@testset "FEM: Nonlinear spring" begin
    @testset "not depending on depth" begin
        pile = FEPileModel(0, 20, 500)
        pile.E .= 1.0e7
        pile.I .= 1.0e-3
        pile.D .= 0.5
        pile.Fext[1] = 1000.0
        k = 2000
        pile.pycurves .= pycurve(y, z) = z > 19 ? 0 : k*y^0.5
        solve!(pile)
        @test log10(pile.u[1]) ≈ 0.1536   atol = 1e-2
        @test log10(pile.u[26]) ≈ -0.0529 atol = 1e-2
        @test log10(abs(pile.θ[1])) ≈ -0.2554  atol = 1e-2
        @test log10(abs(pile.θ[26])) ≈ -0.2964 atol = 1e-2
        @test log10(maximum(abs, pile.M)) ≈ 3.2055 atol = 1e-2
    end
    @testset "depending on depth" begin
        pile = FEPileModel(0, 20, 100)
        pile.E .= 1.0e7
        pile.I .= 1.0e-3
        pile.D .= 0.5
        pile.Fext[1] = 1000
        k = 2000
        pile.pycurves .= pycurve(y, z) = z > 19 ? 0 : k*(19-z)*y^0.5
        solve!(pile)
        @test log10(pile.u[1]) ≈ 0.2391  atol = 1e-3
        @test log10(pile.u[6]) ≈ 0.0385  atol = 1e-3
        @test log10(abs(pile.θ[1])) ≈ -0.1817 atol = 1e-3
        @test log10(abs(pile.θ[6])) ≈ -0.2161 atol = 1e-3
        @test log10(maximum(abs, pile.M)) ≈ 3.3128  atol = 1e-3
    end
end

@testset "solve_disp_load" begin
    # compare with Chang's equation
    pile = FEPileModel(0, 19, 400)
    D = 0.6
    E = 2e8
    I = 0.0002507
    k_h = 3750
    H = 10
    pile.D .= D
    pile.E .= E
    pile.I .= I
    pile.pycurves .= pycurve(y, z) = k_h * y;
    pile.Fext[1] = H
    disp, load = solve_disp_load(pile)
    for (d, f) in zip(disp, load)
        @test d ≈ deflection_chang(f, D, E, I, k_h, 0, 0)  atol = 1e-6
    end
end
