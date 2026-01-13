"""
    ComplexNormal(μ, σ²)

Circularly symmetric complex normal distribution with mean μ and variance σ².
The real and imaginary parts are independent N(Re(μ), σ²/2) and N(Im(μ), σ²/2).
"""
struct ComplexNormal{T<:Real} <: ContinuousUnivariateDistribution
    μ::Complex{T}     # mean
    σ::T              # standard deviation (total, not per component)
    
    function ComplexNormal{T}(μ::Complex{T}, σ::T) where T<:Real
        σ > 0 || throw(ArgumentError("standard deviation must be positive"))
        new{T}(μ, σ)
    end
end

# Convenience constructors
ComplexNormal(μ::Complex{T}, σ::T) where T<:Real = ComplexNormal{T}(μ, σ)
ComplexNormal(μ::T, σ::T) where T<:Real = ComplexNormal{T}(Complex{T}(μ), σ)
ComplexNormal(μ::Number, σ::Real) = ComplexNormal(Complex{Float64}(μ), Float64(σ))
ComplexNormal() = ComplexNormal(0.0 + 0.0im, 1.0)  # Standard complex normal

# Basic properties
Distributions.params(d::ComplexNormal) = (d.μ, d.σ)
Base.eltype(::Type{ComplexNormal{T}}) where T = Complex{T}
Distributions.mean(d::ComplexNormal) = d.μ
Distributions.var(d::ComplexNormal) = abs2(d.σ)


xval(d::ComplexNormal, z::Complex) = muladd(d.σ, z, d.μ)

# Sampling
rand(rng::AbstractRNG, d::ComplexNormal{T}) where T<:Real = xval(d, randn(rng, complex(T)))
function Distributions.rand!(rng::X, d::ComplexNormal{T}, out::AbstractVector{Complex{T}}) where {T<:Real, X<:AbstractRNG}
    randn!(rng, out)
    map!(Base.Fix1(xval, d), out, out)
    return out
end
function Distributions._rand!(rng::X, d::ComplexNormal{T}, out::AbstractVector{Complex{T}}) where {T<:Real, X<:AbstractRNG}
    randn!(rng, out)
    map!(Base.Fix1(xval, d), out, out)
    return out
end

# Log probability density function
function Distributions.logpdf(d::ComplexNormal{T}, z::Number) where T<:Real
    # For circularly symmetric complex normal: p(z) = (1/(πσ²)) * exp(-|z-μ|²/σ²)
    # log p(z) = -log(π) - log(σ²) - |z-μ|²/σ²
    δ = z - d.μ
    return -log(T(π) * abs2(d.σ)) - abs2(δ) / Distributions.var(d)
end

# Probability density function
Distributions.pdf(d::ComplexNormal, z::Number) = exp(Distributions.logpdf(d, z))
Distributions.sampler(d::ComplexNormal{T}) where {T} = d
