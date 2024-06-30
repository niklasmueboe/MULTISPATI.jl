# Reference API

Documentation for `MULTISPATI.jl`'s public interface.

## Index

```@index
Pages = ["api.md"]
```

## API

### MULTISPATI
```@docs
AbstractMultispati
Multispati
fit(::Type{Multispati},
    X::AbstractMatrix{T},
    W::AbstractMatrix{U},
    Q::AbstractMatrix{T}=I,
    D::AbstractMatrix{T}=I / size(X, 2);
    kwargs
) where {T<:Real,U<:Real}
predict(::Multispati, x::AbstractVecOrMat{T}) where {T<:Real}
reconstruct(::Multispati, y::AbstractVecOrMat{T}) where {T<:Real}
moransIbounds
size(::AbstractMultispati)
projection(::AbstractMultispati)
eigvecs(::AbstractMultispati)
eigvals(::AbstractMultispati)
```

### spatialPCA

```@docs
SpatialPCA
fit(::Type{SpatialPCA},
    X::AbstractMatrix{T},
    W::AbstractMatrix{U};
    kwargs
) where {T<:Real,U<:Real}
predict(::SpatialPCA, x::AbstractVecOrMat{T}) where {T<:Real}
reconstruct(::SpatialPCA, y::AbstractVecOrMat{T}) where {T<:Real}
varianceMoransIdecomposition
mean(::SpatialPCA)
```