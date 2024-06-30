module MULTISPATI

import Base: size, show
import LinearAlgebra: eigvals, eigvecs, issymmetric
import MultivariateStats: LinearDimensionalityReduction, projection, reconstruct
import StatsAPI: fit, predict
import Statistics: mean

using Arpack: eigs
using Base: @__dot__
using LinearAlgebra: I, Symmetric, eigen, normalize!
using SparseArrays: issparse

export AbstractMultispati,
    Multispati,
    SpatialPCA,
    fit,
    predict,
    mean,
    eigvals,
    eigvecs,
    projection,
    reconstruct,
    moransIbounds,
    varianceMoransIdecomposition

include("multispati_.jl")
include("spatialPCA.jl")

end # module MULTISPATI
