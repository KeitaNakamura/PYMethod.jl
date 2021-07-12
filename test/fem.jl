function deflection_chang(H::Real, D::Real, E::Real, I::Real, k_h::Real, h::Real, x::Real)
    β = (k_h*D/(4E*I))^(1/4)
    if x < 0
        H / (6E*I*β^3) * (β^3*x^3 + 3β^3*h*x^2 - 3(1+2β*h)*β*x + 3(1+β*h))
    else
        H / (2E*I*β^3) * exp(-β*x) * ((1+β*h)*cos(β*x) - β*h*sin(β*x))
    end
end

@testset "Compare with Chang's equation" begin
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
