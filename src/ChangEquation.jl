struct ChangEquation
    F_t::Float64   # lateral force at top
    M_t::Float64   # bending moment at top
    D::Float64     # diameter of pile
    E::Float64     # Young's modulus of pile
    I::Float64     # second moment of area of pile
    k::Float64     # reaction modulus of the ground
    h::Float64     # height of pile head from the ground
    z_b::Float64   # height of bottom of pile
    z_t::Float64   # height of top of pile
    z_0::Float64   # height of ground surface
    β::Float64
    h_0::Float64
    l_m::Float64   # depth at which M = M_max
    l_0::Float64   # depth of the first fixed point in the ground
    L::Float64     # depth at which deflection angle is zero
    y_t::Float64   # deflection at top
    y_0::Float64   # deflection at ground surface
    θ_t::Float64   # deflection angle at top
    M_max::Float64 # maximum moment in the ground
    function ChangEquation(F_t::Real, M_t::Real, D::Real, E::Real, I::Real, k::Real, h::Real, z_b::Real, z_t::Real, z_0::Real)
        β = sqrt(sqrt(D*k / (4E*I)))
        h_0 = h + M_t / F_t
        l_m = atan(1 / (1+2β*h_0)) / β
        l_0 = atan((1+β*h_0) / (β*h_0)) / β
        L = atan(-(1+2β*h_0)) / β
        y_t = ((1+β*h)^3 + 1/2) / (3E*I*β^3) * F_t + (1+β*h)^2 / (2E*I*β^2) * M_t
        y_0 = (1+β*h_0) / (2E*I*β^3) * F_t
        θ_t = (1+β*h)^2 / (2E*I*β^2) * F_t + (1+β*h) / (E*I*β) * M_t
        M_max = -F_t/2β*sqrt((1+2β*h_0)^2 + 1) * exp(-β*l_m)
        new(F_t, M_t, D, E, I, k, h, z_b, z_t, z_0, β, h_0, l_m, l_0, L, y_t, y_0, θ_t, M_max)
    end
end

function ChangEquation(z_b::Real, z_t::Real; z_0::Real = 0, F_t::Real, M_t::Real = 0, D::Real, E::Real, I::Real, k::Real, head_free::Bool = true)
    h = z_t - z_0
    if !head_free
        β = sqrt(sqrt(D*k / (4E*I)))
        ChangEquation(F_t, -(1+β*h)/2β*F_t, D, E, I, k, h, z_b, z_t, z_0)
    else
        ChangEquation(F_t, M_t, D, E, I, k, h, z_b, z_t, z_0)
    end
end

function calculate_deflection(eq::ChangEquation, z::Real)
    @assert eq.z_b ≤ z ≤ eq.z_t
    E = eq.E
    I = eq.I
    β = eq.β
    if z > eq.z_0
        x = -(z - eq.z_t)
        eq.y_t - eq.θ_t*x + eq.M_t/(2E*I)*x^2 + eq.F_t/(6E*I)*x^3
    else
        x = -(z - eq.z_0)
        h_0 = eq.h_0
        eq.F_t/(2E*I*β^3) * exp(-β*x) * ((1+β*h_0) * cos(β*x) - β*h_0*sin(β*x))
    end
end

function calculate_moment(eq::ChangEquation, z::Real)
    @assert eq.z_b ≤ z ≤ eq.z_t
    x = -(z - eq.z_0)
    β = eq.β
    if z > eq.z_0
        x = -(z - eq.z_t)
        -eq.M_t - eq.F_t*x
    else
        x = -(z - eq.z_0)
        h_0 = eq.h_0
        -eq.F_t/β * exp(-β*x) * (β*h_0*cos(β*x) + (1+β*h_0)*sin(β*x))
    end
end

function calculate_shearforce(eq::ChangEquation, z::Real)
    @assert eq.z_b ≤ z ≤ eq.z_t
    x = -(z - eq.z_0)
    β = eq.β
    if z > eq.z_0
        x = -(z - eq.z_t)
        -eq.F_t
    else
        x = -(z - eq.z_0)
        h_0 = eq.h_0
        -eq.F_t * exp(-β*x) * (cos(β*x) - (1+2β*h_0)*sin(β*x))
    end
end
