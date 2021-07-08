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
function reset!(x::FEMVector{T}) where {T}
    fill!(x.data, zero(T))
    fill!(x.mask, false)
    x
end

struct Beam{T}
    U::FEMVector{T}
    inds::UnitRange{Int}
    # constants
    l::T
    E::T
    I::T
end

struct FEPileModel{T} <: AbstractVector{Beam{T}}
    coordinates::LinRange{T}
    # global vectors
    U::FEMVector{T}
    Fext::Vector{T}
    # parameters
    E::Vector{T}
    I::Vector{T}
    D::Vector{T}
    # p-y curves
    pycurves::Vector{Any}
end

function FEPileModel(bottom::Real, top::Real, nelements::Int)
    coords = LinRange(top, bottom, nelements + 1)
    n = length(coords)
    ndofs = n * 2 # one direction, TODO: handle two directions
    T = eltype(coords)
    # global vectors and matrix
    U = FEMVector{T}(ndofs)
    Fext = zeros(T, ndofs)
    # parameters
    Eᵢ = zeros(T, n)
    Iᵢ = zeros(T, n)
    Dᵢ = zeros(T, n)
    # p-y curves
    pycurves = Vector{Any}(undef, n)
    FEPileModel(coords, U, Fext, Eᵢ, Iᵢ, Dᵢ, pycurves)
end

function Base.getproperty(model::FEPileModel, name::Symbol)
    U = getfield(model, :U)
    Fext = getfield(model, :Fext)
    if name == :u
        view(U, 1:2:length(U))
    elseif name == :θ
        view(U, 2:2:length(U))
    elseif name == :f
        view(Fext, 1:2:length(Fext))
    elseif name == :M
        view(Fext, 2:2:length(Fext))
    else
        getfield(model, name)
    end
end

@generated Base.propertynames(model::FEPileModel) =
    :(($(map(QuoteNode, fieldnames(model))...), :u, :θ, :F, :M))

Base.size(model::FEPileModel) = (length(model.coordinates)-1,)
function Base.getindex(model::FEPileModel, i::Int)
    @boundscheck checkbounds(model, i)
    Z = model.coordinates
    l = Z[i+1] - Z[i]
    E = mean(model.E[i:i+1])
    I = mean(model.I[i:i+1])
    ind = 2(i-1) + 1
    Beam(model.U, ind:ind+3, abs(l), E, I)
end

function stiffness_matrix(l::Real, E::Real, I::Real)
    E*I/l^3 * [ 12  6l   -12  6l
                6l  4l^2 -6l  2l^2
               -12 -6l    12 -6l
                6l  2l^2 -6l  4l^2]
end
stiffness_matrix(beam::Beam) = stiffness_matrix(beam.l, beam.E, beam.I)

function reaction_forces(beam::Beam)
    inds = beam.inds
    Fext = stiffness_matrix(beam) * beam.U[inds]
    f = @view Fext[1:2:end]
    M = @view Fext[2:2:end]
    (; f, M)
end

function distributions(model::FEPileModel{T}) where {T}
    F = T[]
    M = T[]
    for i in 1:length(model)
        beam = model[i]
        F_beam, M_beam = reaction_forces(beam)
        push!(F, F_beam[1])
        push!(M, M_beam[1])
        if i == length(model)
            push!(F, -F_beam[2])
            push!(M, -M_beam[2])
        end
    end
    (; model.u, model.θ, F, M)
end

function stiffness_pycurve(pycurve, D::Real, l′::Real, y::Real, z::Real)
    k = ForwardDiff.derivative(y -> oftype(y, pycurve(y, z)), y)
    k * D * l′
end

function soil_reaction_force(pycurve, D::Real, l′::Real, y::Real, z::Real)
    pycurve(y, z) * D * l′ # (pressure) * (area)
end

function assemble_force_vector!(Fint::AbstractVector, U::AbstractVector, model::FEPileModel{T}) where {T}
    Z = model.coordinates
    for beam in model
        inds = beam.inds
        Fint[inds] += stiffness_matrix(beam) * U[inds]
    end
    for i in 1:length(Z)
        ind = 2(i-1) + 1
        y = U[i]
        z = Z[i]
        D = model.D[i]
        if i == 1
            l′ = mean(Z[i:i+1])
        elseif i == length(Z)
            l′ = mean(Z[i-1:i])
        else
            l′ = mean(Z[i:i+1]) - mean(Z[i-1:i])
        end
        Fint[ind] += soil_reaction_force(model.pycurves[i], D, l′, y, z)
    end
end

function solve!(model::FEPileModel{T}) where {T}
    U = parent(model.U)
    Fext = model.Fext
    Fint = similar(Fext)
    K = similar(Fint, length(Fext), length(Fext))
    fdofs = .!model.U.mask
    fdofs[end-1:end] .= false # bottom fixed
    for i in 1:20
        fill!(Fint, zero(T))
        ForwardDiff.jacobian!(K, (y, x) -> assemble_force_vector!(y, x, model), Fint, U)
        R = Fint - Fext
        norm(R[fdofs]) < 1e-8 && return
        U[fdofs] .-= K[fdofs, fdofs] \ R[fdofs]
    end
    error("too mach iterations")
end

function reset!(model::FEPileModel{T}) where {T}
    reset!(model.U)
    fill!(model.Fext, zero(T))
    model
end

@recipe function f(model::FEPileModel)
    xguide --> "Lateral displacement"
    yguide --> "Coordinates"
    label --> ""
    (model.u, model.coordinates)
end
