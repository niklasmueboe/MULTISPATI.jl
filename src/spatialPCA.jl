
"""
    SpatialPCA <: AbstractMultispati
"""
struct SpatialPCA{T<:Real} <: AbstractMultispati
    mean
    proj::AbstractMatrix{T}
    eigvals::AbstractVector{T}
    W::AbstractMatrix{T}
end

"""
    mean(M::SpatialPCA)

Returns the mean vector (of length `d`).
"""
mean(M::SpatialPCA) = M.mean

"""
    predict(M::SpatialPCA, x::AbstractVecOrMat{<:Real})

Transform the observations `x` with the SpatialPCA model `M`.

Here, `x` can be either a vector of length `d` or a matrix where each column is an observation.
"""
predict(M::SpatialPCA, x::AbstractVecOrMat{T}) where {T<:Real} =
    transpose(projection(M)) * centralize(M, x)

"""
    reconstruct(M::SpatialPCA, y::AbstractVecOrMat{<:Real})

Approximately reconstruct the observations `y` to the original space using the SpatialPCA model `M`.

Here, `y` can be either a vector of length `p` or a matrix where each column
gives the components for an observation.
"""
reconstruct(M::SpatialPCA, y::AbstractVecOrMat{T}) where {T<:Real} =
    decentralize(M, projection(M) * y)

function centralize(M::SpatialPCA, x)
    if isnothing(mean(M)) || issparse(x)
        return x
    else
        return x .- mean(M)
    end
end

function decentralize(M::SpatialPCA, x)
    if isnothing(mean(M)) || issparse(x)
        return x
    else
        return x .+ mean(M)
    end
end

"""
    varianceMoransIdecomposition(M::SpatialPCA, X)

Decompose the eigenvalues into a variance and Moran's I contribution given the model `M` 
and matrix `X` which was used for fitting the model.
"""
function varianceMoransIdecomposition(M::SpatialPCA, X)
    transformedX = predict(M, X)
    laggedX = transformedX * M.W
    w = 1 / size(M.W, 1) # sum of row_weights but because its normed the sum is n
    variance = sum(transformedX .* transformedX .* w; dims=2)
    moransI = sum(transformedX .* laggedX; dims=2) ./ variance # TODO ? 
    return variance, moransI
end

function show(io::IO, M::SpatialPCA)
    idim, odim = size(M)
    return print(io, "SpatialPCA(indim = $idim, outdim = $odim)")
end

"""
    fit(SpatialPCA, X, W; ...)

Perform spatialPCA over the data given a matrix `X` and `W`. Each column of `X` is an **observation**. 
`W` is a connectivity matrix where ``w_{ij}`` is the connection from j -> i.

**Keyword arguments**

- `maxoutdim`: The output dimension, i.e. dimension of the transformed space (*min(d, nc-1)*)
- `solver`: The choice of solver:
    - `:eig`: uses `LinearAlgebra.eigen` (*default*)
    - `:eigs`: uses `Arpack.eigs` (always used for sparse data)
- `tol`: Convergence tolerance for `eigs` solver (*default* `0.0`)
- `maxiter`: Maximum number of iterations for `eigs` solver (*default* `300`)
- `center_sparse`: Center sparse matrix `X` (dense X will always be centered) (*default* `false`)

**References**

[T. Jombart, et al. "Revealing cryptic spatial patterns in genetic variability by a new multivariate method."
*Heredity* (2008)](https://doi.org/10.1038/hdy.2008.34)
"""
function fit(
    ::Type{SpatialPCA},
    X::AbstractMatrix{T},
    W::AbstractMatrix{U};
    maxoutdim::Union{Nothing,Integer,Tuple{Integer,Integer}}=nothing,
    solver::Symbol=:eig,
    tol::Real=0.0,
    maxiter::Real=300,
    center_sparse::Bool=false,
) where {T<:Real,U<:Real}
    d, n = size(X)
    m1, m2 = size(W)

    if m1 != m2
        throw(DimensionMismatch("`W` must be a square Matrix"))
    end
    if m1 != n
        throw(DimensionMismatch("# cols of `X` must match dimensions of `W`"))
    end

    if center_sparse && issparse(X)
        X = Matrix(X)
    end

    if !issparse(X)
        meanvec = vec(mean(X; dims=2))
        X = X .- meanvec
    else
        meanvec = nothing
    end

    normalize!.(eachcol(W), 1)

    maxoutdim = validate_maxoutdim(d, n, maxoutdim)

    H = (X * (transpose(W) + W) * transpose(X)) / 2n

    @assert isapprox(H, H') # issymmetric does not have option to specify tolerance

    eigenvals, eigenvecs = multispati_decomposition(
        Symmetric(H), maxoutdim; solver=solver, tol=tol, maxiter=maxiter
    )

    return SpatialPCA(mean, eigenvecs, eigenvals, W)
end
