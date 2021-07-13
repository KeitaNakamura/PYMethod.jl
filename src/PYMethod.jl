module PYMethod

using LinearAlgebra, Statistics
using StaticArrays, ForwardDiff
using RecipesBase

using Base: @_propagate_inbounds_meta

export
    ChangEquation,
    FEPileModel,
    solve!

include("utils.jl")
include("ChangEquation.jl")
include("fem.jl")

end # module
