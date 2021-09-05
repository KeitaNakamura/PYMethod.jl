module PYMethod

using LinearAlgebra, Statistics, SparseArrays
using StaticArrays, ForwardDiff
using RecipesBase

using Base: @_propagate_inbounds_meta

export
    ChangEquation,
    FEPileModel,
    solve!

fillzero!(x) = fill!(x, zero(eltype(x)))

include("ChangEquation.jl")
include("fem.jl")

end # module
