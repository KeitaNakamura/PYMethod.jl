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

struct FEPileModel{T}
    coordinates::LinRange{T}
    # global vectors and matrix (Ax = b)
    x::FEMVector{T}
    b::Vector{T}
    A::Matrix{T}
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
    x = FEMVector{T}(ndofs)
    b = zeros(T, ndofs)
    A = zeros(T, ndofs, ndofs)
    # parameters
    Eᵢ = zeros(T, n)
    Iᵢ = zeros(T, n)
    Dᵢ = zeros(T, n)
    # p-y curves
    pycurves = Vector{Any}(undef, n)
    FEPileModel(coords, x, b, A, Eᵢ, Iᵢ, Dᵢ, pycurves)
end

function Base.getproperty(model::FEPileModel, name::Symbol)
    x = getfield(model, :x)
    b = getfield(model, :b)
    if name == :u
        view(x, 1:2:length(x))
    elseif name == :θ
        view(x, 2:2:length(x))
    elseif name == :F
        view(b, 1:2:length(b))
    elseif name == :M
        view(b, 2:2:length(b))
    else
        getfield(model, name)
    end
end

@generated Base.propertynames(model::FEPileModel) =
    :(($(map(QuoteNode, fieldnames(model))...), :u, :θ, :F, :M))

function stiffness_matrix(l::Real, E::Real, I::Real)
    E*I/l^3 * [ 12  6l   -12  6l
                6l  4l^2 -6l  2l^2
               -12 -6l    12 -6l
                6l  2l^2 -6l  4l^2]
end

function stiffness_pycurve(pycurve, D::Real, l′::Real, y::Real, z::Real)
    k = ForwardDiff.derivative(y -> oftype(y, pycurve(y, z)), y)
    k * D * l′
end

function assemble_stiffness_matrix!(model::FEPileModel{T}) where {T}
    K = fill!(model.A, zero(T))
    Z = model.coordinates
    for i in 1:length(Z)-1
        ind = 2(i-1) + 1
        l = Z[i+1] - Z[i]
        E = mean(model.E[i:i+1])
        I = mean(model.I[i:i+1])
        K[ind:ind+3, ind:ind+3] += stiffness_matrix(l, E, I)
    end
    for i in 1:length(Z)-1
        ind = 2(i-1) + 1
        u = model.u[i]
        z = Z[i]
        D = mean(model.D[i:i+1])
        if i == 1
            l′ = mean(Z[i:i+1])
        elseif i == length(Z)
            l′ = mean(Z[i-1:i])
        else
            l′ = mean(Z[i:i+1]) - mean(Z[i-1:i])
        end
        K[ind, ind] += stiffness_pycurve(model.pycurves[i], D, l′, u, z)
    end
end

function solve!(u::AbstractVector, F::AbstractVector, K::AbstractMatrix, fdofs::BitVector)
    # sweep
    @. fdofs = !fdofs
    F .-= K[:, fdofs] * u[fdofs]
    # solve equation
    @. fdofs = !fdofs
    u[fdofs] = K[fdofs, fdofs] \ F[fdofs]
    u
end

function solve!(model::FEPileModel)
    u = parent(model.x)
    u_n = similar(u)
    F = similar(model.b)
    fdofs = .!model.x.mask
    fdofs[end-1:end] .= false # bottom fixed
    while true
        copy!(u_n, u)
        copy!(F, model.b)
        assemble_stiffness_matrix!(model)
        solve!(u, F, model.A, fdofs)
        norm(u - u_n) < 1e-8 && break
    end
    model.u
end

function reset!(model::FEPileModel{T}) where {T}
    reset!(model.x)
    fill!(model.b, zero(T))
    model
end
