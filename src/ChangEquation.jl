"""
    ChangEquation(top, bottom; parameters...)

# Parameters
* `z_0`: height of ground surface (`z_0 = 0` by default)
* `F_t`: Lateral load at pile head
* `M_t`: Moment at pile head (`M_t = 0` by default)
* `D`: Diameter of pile
* `E`: Young's modulus of pile
* `I`: Second moment of area of pile
* `k`: Modulus of subgrade reaction

# Examples
```jldoctest
julia> eq = ChangEquation(-19, 1; F_t = 10, D = 0.6, E = 2e8, I = 0.0002507, k = 3750);
```
"""
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

function ChangEquation(z_t::Real, z_b::Real; z_0::Real = 0, F_t::Real, M_t::Real = 0, D::Real, E::Real, I::Real, k::Real, head_free::Bool = true)
    h = z_0 - z_t
    if !head_free
        β = sqrt(sqrt(D*k / (4E*I)))
        ChangEquation(F_t, -(1+β*h)/2β*F_t, D, E, I, k, h, z_b, z_t, z_0)
    else
        ChangEquation(F_t, M_t, D, E, I, k, h, z_b, z_t, z_0)
    end
end

"""
    calculate_deflection(::ChangEquation, z)

Calculate deflection at depth `z`.
"""
function calculate_deflection(eq::ChangEquation, z::Real)
    @assert eq.z_t ≤ z ≤ eq.z_b
    E = eq.E
    I = eq.I
    β = eq.β
    if z > eq.z_0
        x = z - eq.z_0
        h_0 = eq.h_0
        eq.F_t/(2E*I*β^3) * exp(-β*x) * ((1+β*h_0) * cos(β*x) - β*h_0*sin(β*x))
    else
        x = z - eq.z_t
        eq.y_t - eq.θ_t*x + eq.M_t/(2E*I)*x^2 + eq.F_t/(6E*I)*x^3
    end
end

"""
    calculate_moment(::ChangEquation, z)

Calculate moment at height `z`.
"""
function calculate_moment(eq::ChangEquation, z::Real)
    @assert eq.z_t ≤ z ≤ eq.z_b
    β = eq.β
    if z > eq.z_0
        x = z - eq.z_0
        h_0 = eq.h_0
        -eq.F_t/β * exp(-β*x) * (β*h_0*cos(β*x) + (1+β*h_0)*sin(β*x))
    else
        x = z - eq.z_t
        -eq.M_t - eq.F_t*x
    end
end

"""
    calculate_shearforce(::ChangEquation, z)

Calculate shear force at height `z`.
"""
function calculate_shearforce(eq::ChangEquation, z::Real)
    @assert eq.z_t ≤ z ≤ eq.z_b
    β = eq.β
    if z > eq.z_0
        x = z - eq.z_0
        h_0 = eq.h_0
        -eq.F_t * exp(-β*x) * (cos(β*x) - (1+2β*h_0)*sin(β*x))
    else
        x = z - eq.z_t
        -eq.F_t
    end
end

@recipe function f(eq::ChangEquation)
    label --> ""
    layout := (1, 3)
    coords = LinRange(eq.z_t, eq.z_b, 1000)
    u = calculate_deflection.(Ref(eq), coords)
    M = calculate_moment.(Ref(eq), coords)
    F = calculate_shearforce.(Ref(eq), coords)
    @series begin
        subplot := 1
        xguide --> "Deflection"
        yguide --> "Coordinate"
        (u, coords)
    end
    @series begin
        subplot := 2
        xguide --> "Moment"
        (M, coords)
    end
    @series begin
        subplot := 3
        xguide --> "Shear force"
        (F, coords)
    end
end
