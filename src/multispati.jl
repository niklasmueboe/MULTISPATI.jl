# MULTISPATI

# TODO: :auto solver
# TODO: centering of Multispati?

"""
    AbstractMultispati
"""
abstract type AbstractMultispati <: LinearDimensionalityReduction end

"""
    Multispati <: AbstractMultispati
"""
struct Multispati{T<:Real} <: AbstractMultispati
    proj::AbstractMatrix{T}
    eigvals::AbstractVector{T}
    W::AbstractMatrix{T}
end

"""
    size(M)

Returns a tuple with the dimensions of input (the dimension of the observation space)
and output (the dimension of the principal subspace).
"""
size(M::AbstractMultispati) = size(M.proj)

"""
    projection(M::AbstractMultispati)

Returns the projection matrix (of size `(d, p)`). Each column of the projection matrix corresponds to a eigenvector.
The eigenvectors are arranged in ascending order of the eigenvalues.
"""
projection(M::AbstractMultispati) = M.proj

"""
    eigvecs(M::AbstractMultispati)

Get the eigenvectors of the Multispati model `M`.
"""
eigvecs(M::AbstractMultispati) = projection(M)

"""
    eigvals(M::AbstractMultispati)

Get the eigenvalues of the Multispati model `M`.
"""
eigvals(M::AbstractMultispati) = M.eigvals

"""
    predict(M::Multispati, x::AbstractVecOrMat{<:Real})

Transform the observations `x` with the Multispati model `M`.

Here, `x` can be either a vector of length `d` or a matrix where each column is an observation.
"""
predict(M::Multispati, x::AbstractVecOrMat{T}) where {T<:Real} =
    transpose(projection(M)) * x

"""
    reconstruct(M::Multispati, y::AbstractVecOrMat{<:Real})

Approximately reconstruct the observations `y` to the original space using the Multispati model `M`.

Here, `y` can be either a vector of length `p` or a matrix where each column
gives the components for an observation.
"""
reconstruct(M::Multispati, y::AbstractVecOrMat{T}) where {T<:Real} = projection(M) * y

"""
    moransIbounds(M::AbstractMultispati; sparse_approx::Bool=true)

Return the bounds and expected value for Moran's I given the model `M` in the order 
``I_{min}``, ``I_{max}``, ``I_0`` 
"""
function moransIbounds(M::AbstractMultispati; sparse_approx::Bool=true)
    function doublecenter!(X)
        rowmeans = mean(W; dims=1)
        colmeans = mean(W; dims=2) - mean(rowmeans)
        @. X -= rowmeans - colmeans
    end

    L = (M.W + transpose(M.W)) / 2
    n = size(L, 1)

    if !sparse_approx && issparse(L)
        L = Matrix(L)
    end

    if !issparse(L)
        doublecenter!(L)
    end

    eigenvals = n / sum(L) * eigs(L; nev=2, which=:BE, ritzvec=false)

    I0 = -1 / (n - 1)
    Imin = min(eigenvals)
    Imax = max(eigenvals)

    return Imin, Imax, I0
end

function show(io::IO, M::Multispati)
    idim, odim = size(M)
    return print(io, "Multispati(indim = $idim, outdim = $odim)")
end

"""
    fit(Multispati, X, W, Q=I, D=I / size(X, 2); ...)

Perform Multispati over the data given a matrix `X`. Each column of `X` is an **observation**.
`W` is a connectivity matrix where ``w_{ij}`` is the connection from j -> i.
`Q` is a symmetric matrix of size `n` and `D` a symmetric matrix of size `d`

**Keyword arguments**

- `maxoutdim`: The output dimension, i.e. dimension of the transformed space (*min(d, nc-1)*)
- `solver`: The choice of solver:
    - `:eig`: uses `LinearAlgebra.eigen` (*default*)
    - `:eigs`: uses `Arpack.eigs` (always used for sparse data)
- `tol`: Convergence tolerance for `eigs` solver (*default* `0.0`)
- `maxiter`: Maximum number of iterations for `eigs` solver (*default* `300`)

**References**

[S. Dray, et al. "Spatial ordination of vegetation data using a generalization of Wartenberg's 
multivariate spatial correlation." *Journal of vegetation science* (2008)](https://doi.org/10.3170/2007-8-18312)

[de la Cruz and Holmes. "The duality diagram in data analysis: Examples of modern applications." 
*The annals of applied statistics* (2011)](https://doi.org/10.1214/10-aoas408)
"""
function fit(
    ::Type{Multispati},
    X::AbstractMatrix{T},
    W::AbstractMatrix{U},
    Q::AbstractMatrix{T}=I,
    D::AbstractMatrix{T}=I / size(X, 2);
    maxoutdim::Union{Nothing,Integer,Tuple{Integer,Integer}}=nothing,
    solver::Symbol=:eig,
    tol::Real=0.0,
    maxiter::Real=300,
) where {T<:Real,U<:Real}
    m, n = size(X)
    w1, w2 = size(W)

    if w1 != w2
        throw(DimensionMismatch("`W` must be a square Matrix"))
    end

    if !issymmetic(Q)
        throw(ArgumentError("`Q` must be symmetric"))
    end
    if !issymmetic(D)
        throw(ArgumentError("`D` must be symmetric"))
    end

    if w1 != n
        throw(DimensionMismatch("# cols of `X` must match dimensions of `W`"))
    end
    if d1 != n
        throw(DimensionMismatch("# cols of `X` must match dimensions of `D`"))
    end
    if q1 != m
        throw(DimensionMismatch("# rows of `X` must match dimensions of `Q`"))
    end

    normalize!.(eachcol(W), p=1)

    H = (X * (W * D + D * transpose(W)) * transpose(X) * Q) / 2

    eigenvals, eigenvecs = multispati_decomposition(
        H, maxoutdim; solver=solver, tol=tol, maxiter=maxiter
    )

    return Multispati(eigenvecs, eigenvals, W)
end

function validate_maxoutdim(d::Integer, n::Integer, maxoutdim)
    maxout = max(d, n)
    if !isnothing(maxoutdim)
        if maxoutdim isa Integer && maxoutdim >= maxout
            return nothing
        elseif maxoutdim isa Tuple && sum(maxoutdim) >= maxout
            return nothing
        end
    end
    return maxoutdim
end

# TODO: return negative instead of lowest eigenvalues
function multispati_decomposition(H, outdim; kwargs...)
    outdim = validate_maxoutdim(m, n, outdim)

    # TODO: remove zero eigenvalues?
    if solver == :eigs || issparse(H)
        eigenvals, eigenvecs = arpack_decomposition(H, outdim...; kwargs...)
        return real.(eigenvals), real.(eigenvecs)
    elseif solver == :eig
        return la_decomposition(Symmetric(H), outdim...)
    else
        throw(ArgumentError("Invalid solver name $(solver)"))
    end
end

# ARPACK solver
function arpack_decomposition(H, ::Nothing; kwargs...)
    return eigs(H; nev=size(H, 1) - 1, which=:LM, ritzvec=true, kwargs...)
end
function arpack_decomposition(H, npos; kwargs...)
    return eigs(H; nev=npos, which=:LR, ritzvec=true, kwargs...)
end
function arpack_decomposition(H, npos, nneg; kwargs...)
    if npos == 0
        return eigs(H; nev=nneg, which=:SR, ritzvec=true, kwargs...)
    elseif nneg == 0
        return eigs(H; nev=npos, which=:LR, ritzvec=true, kwargs...)
    else
        #TODO
        return eigs(H; nev=npos, which=:BE, ritzvec=true, kwargs...)
    end
end

# LinearAlgebra solver
function la_decomposition(H, npos::Union{Nothing,Integer})
    val, vec = eig(H)
    n = size(val)
    if !isnothing(npos)
        ind = (n - npos + 1):n
        val = val[ind]
        vec = vec[:, ind]
    end
    return val, vec
end
function la_decomposition(H, npos::Integer, nneg::Integer)
    val, vec = eig(H)
    n = size(val)
    if nneg == 0
        ind = (n - npos + 1):n
    elseif npos == 0
        ind = 1:nneg
    else
        ind = (1:nneg; (n - npos + 1):n)
    end

    val = val[ind]
    vec = vec[:, ind]
    return val, vec
end
