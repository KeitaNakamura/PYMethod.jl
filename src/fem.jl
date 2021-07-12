struct FEMVector{T} <: AbstractVector{T}
    data::Vector{T}
    mask::BitVector
end

FEMVector{T}(n::Int) where {T} = FEMVector(zeros(T, n), falses(n))
Base.parent(x::FEMVector) = x.data
Base.size(x::FEMVector) = size(x.data)
Base.getindex(x::FEMVector, i::Int) = (@_propagate_inbounds_meta; x.data[i])
function Base.setindex!(x::FEMVector, v, i::Int)
    @_propagate_inbounds_meta
    x.data[i] = v
    x.mask[i] = true
    x
end
function reset!(x::FEMVector)
    fillzero!(x.data)
    fillzero!(x.mask)
    x
end

struct Beam{T}
    inds::SVector{4, Int}
    # constants
    l::T
    E::T
    I::T
end

struct FEPileModel{T} <: AbstractVector{Beam{T}}
    coordinates::LinRange{T}
    # global vectors
    U::FEMVector{T}
    K::Matrix{T}
    Bext::Vector{T}
    # parameters
    E::Vector{T}
    I::Vector{T}
    D::Vector{T}
    # p-y curves
    pycurves::Vector{Any}
end

"""
    FEPileModel(bottom::Real, top::Real, nelements::Int) -> pile

Construct an object of the finite element model to simulate lateral behavior of pile.
The `i`th `Beam` element can be accessed by `pile[i]`.
The following values are vectors storing each nodal value.

# Parameters

* `pile.E`: Young's modulus of pile
* `pile.I`: Second moment of area of pile
* `pile.D`: Diameter of pile

# Boundary conditions

* `pile.u`: Lateral displacement
* `pile.θ`: Angle of deflection (where `θ` can be typed by `\\theta<tab>`)
* `pile.Fext`: External lateral force
* `pile.Mext`: External moment

The above vectors will be updated after `solve!` the problem.

# Nodal variables

* `pile.coordinates`: Coordinates of nodes
* `pile.F`: Internal lateral force
* `pile.M`: Internal moment

The internal lateral force and moment vectors will be updated after `solve!` the problem.
"""
function FEPileModel(bottom::Real, top::Real, nelements::Int)
    coords = LinRange(top, bottom, nelements + 1)
    n = length(coords)
    ndofs = n * 2 # one direction, TODO: handle two directions
    T = eltype(coords)
    # global vectors and matrix
    U = FEMVector{T}(ndofs)
    Bext = zeros(T, ndofs)
    K = zeros(T, ndofs, ndofs)
    # parameters
    Eᵢ = zeros(T, n)
    Iᵢ = zeros(T, n)
    Dᵢ = zeros(T, n)
    # p-y curves
    pycurves = Vector{Any}(undef, n)
    pycurves .= (pycurve(y, z) = zero(T))
    FEPileModel(coords, U, K, Bext, Eᵢ, Iᵢ, Dᵢ, pycurves)
end

function Base.getproperty(model::FEPileModel, name::Symbol)
    U = getfield(model, :U)
    Bext = getfield(model, :Bext)
    if name == :u
        view(U, 1:2:length(U))
    elseif name == :θ
        view(U, 2:2:length(U))
    elseif name == :Fext
        view(Bext, 1:2:length(Bext))
    elseif name == :Mext
        view(Bext, 2:2:length(Bext))
    elseif name == :F
        internal_force(model, 1)
    elseif name == :M
        internal_force(model, 2)
    else
        getfield(model, name)
    end
end

@generated Base.propertynames(model::FEPileModel) =
    :(($(map(QuoteNode, fieldnames(model))...), :u, :θ, :Fext, :Mext, :F, :M))

Base.size(model::FEPileModel) = (length(model.coordinates)-1,)
function Base.getindex(model::FEPileModel, i::Int)
    @boundscheck checkbounds(model, i)
    Z = model.coordinates
    l = Z[i+1] - Z[i]
    E = mean(@view model.E[i:i+1])
    I = mean(@view model.I[i:i+1])
    ind = 2(i-1) + 1
    Beam(SVector{4}(ind:ind+3), abs(l), E, I)
end

function stiffness_matrix(l::Real, E::Real, I::Real)
    E*I/l^3 * @SMatrix [ 12  6l   -12  6l
                         6l  4l^2 -6l  2l^2
                        -12 -6l    12 -6l
                         6l  2l^2 -6l  4l^2]
end
"""
    stiffness_matrix(::Beam)

Construct element stiffness matrix from `Beam`.
"""
stiffness_matrix(beam::Beam) = stiffness_matrix(beam.l, beam.E, beam.I)

function mass_matrix(l::Real)
    @SMatrix [ 13l/35    -11l^2/210   9l/70    13l^2/420
              -11l^2/210    l^3/105 -13l^2/420  -l^3/140
                9l/70    -13l^2/420  13l/35    11l^2/210
               13l^2/420   -l^3/140  11l^2/210   l^3/105]
end
"""
    mass_matrix(::Beam)

Construct mass matrix from `Beam`.
"""
mass_matrix(beam::Beam) = mass_matrix(beam.l)

function internal_force(model::FEPileModel{T}, id::Int) where {T}
    F = Vector{T}(undef, length(model)+1)
    for i in 1:length(model)
        beam = model[i]
        inds = beam.inds
        Bext = stiffness_matrix(beam) * model.U[inds]
        F_beam = Bext[id:2:end]
        F[i] = F_beam[1]
        if i == length(model)
            F[i+1] = -F_beam[2]
        end
    end
    F
end

function assemble_force_vector!(Fint::AbstractVector, U::AbstractVector, model::FEPileModel)
    P = zero(U)
    Z = model.coordinates
    for beam in model
        inds = beam.inds
        Fint[inds] += stiffness_matrix(beam) * U[inds]
    end
    for i in 1:length(Z)
        ind = 2(i-1) + 1
        y = U[ind]
        z = Z[i]
        D = model.D[i]
        P[ind] = D * model.pycurves[i](y, z)
    end
    for beam in model
        inds = beam.inds
        Fint[inds] += mass_matrix(beam) * P[inds]
    end
end

"""
    solve!(::FEPileModel)

Solve the finite element problem.
"""
function solve!(model::FEPileModel{T}) where {T}
    U = parent(model.U)
    K = model.K
    Bext = model.Bext
    Fint = similar(Bext)
    fdofs = .!model.U.mask
    fdofs[end-1:end] .= false # bottom fixed
    residuals = T[]
    for i in 1:20
        fillzero!(Fint)
        ForwardDiff.jacobian!(K, (y, x) -> assemble_force_vector!(y, x, model), Fint, U)
        R = Fint - Bext
        push!(residuals, norm(R[fdofs]))
        residuals[end] < 1e-6 && return (Bext .= Fint; residuals)
        U[fdofs] .-= K[fdofs, fdofs] \ R[fdofs]
    end
    Bext .= Fint
    residuals
end

"""
    clear_boundary_conditions!(::FEPileModel)

Clear boundary conditions.
The parameters and p-y curves are remained.
"""
function clear_boundary_conditions!(model::FEPileModel)
    reset!(model.U)
    fillzero!(model.Bext)
    model
end

@recipe function f(model::FEPileModel)
    label --> ""
    layout := (1, 4)
    u = model.u
    θ = model.θ
    F = model.F
    M = model.M
    @series begin
        subplot := 1
        xguide --> "Lateral displacement"
        yguide --> "Coordinates"
        (u, model.coordinates)
    end
    @series begin
        subplot := 2
        xguide --> "Angle"
        (θ, model.coordinates)
    end
    @series begin
        subplot := 3
        xguide --> "Lateral force"
        (F, model.coordinates)
    end
    @series begin
        subplot := 4
        xguide --> "Moment"
        (M, model.coordinates)
    end
end

#=
julia> pile = FEPileModel(0, 10, 50);
julia> pile.E .= 2e8;
julia> pile.D .= 0.6;
julia> pile.I .= 0.0002507;
julia> pile.pycurves .= pycurve(y, z) = z > 8 ? 0 : 3750*(8-z)*y;
julia> pile.Fext[1] = 10;
julia> solve!(pile);
julia> plot(pile)

julia> res = [ 1.06737867E+00
        5.32945027E-01
        4.85644953E-01
        4.40098942E-01
        3.96458642E-01
        3.54863263E-01
        3.15435127E-01
        2.78276340E-01
        2.43466519E-01
        2.11061439E-01
        1.81092521E-01
        1.53567040E-01
        1.28468952E-01
        1.05760224E-01
        8.53825843E-02
        6.72595680E-02
        5.12987935E-02
        3.73943655E-02
        2.54293357E-02
        1.52781513E-02
        6.80903016E-03
       -1.13786859E-04
       -5.62795038E-03
       -9.87111600E-03
       -1.29792826E-02
       -1.50853724E-02
       -1.63180635E-02
       -1.68008765E-02
       -1.66515144E-02
       -1.59814456E-02
       -1.48957174E-02
       -1.34929830E-02
       -1.18657200E-02
       -1.01006162E-02
       -8.27909513E-03
       -6.47795175E-03
       -4.77006434E-03
       -3.22514903E-03
       -1.91051851E-03
       -8.91805061E-04
       -2.33605204E-04
        0.00000000E+00
       ]
julia> ys = [10, (8:-0.2:0)...]

julia> plot(res*0.01, ys)
julia> plot!(pile.u, pile.coordinates)
=#