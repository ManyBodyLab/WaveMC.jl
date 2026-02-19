""" 
    AbstractWavefunction

This abstract type defines the interface for wavefunction objects used in Monte Carlo simulations.
Concrete implementations must provide the following methods:
- (psi)(x): Evaluates the wavefunction at the given coordinates x.
- inputtype(YourWavefunctionType): Returns the data type of the input coordinates (e.g., Float64, ComplexF64, Tuple{Float64, Float64}).
- inputlength(YourWavefunctionType): Returns the number of coordinates the wavefunction takes as input.

Optional methods are:
- coordinate_projector(YourWavefunctionType): Returns a function that projects coordinates into the valid domain for the wavefunction. Defaults to the identity function.
- coordinate_transformer(YourWavefunctionType): Returns a function that transforms coordinates before they are used in the wavefunction. Defaults to the identity function.

Methods which should be overloaded for improved performance:
- logdensity(psi::YourWavefunctionType, x): Computes the logarithm of the probability density at coordinates x, corresponds to log(abs2(psi(x))).
- delta_logdensity(psi::YourWavefunctionType, xnew::AbstractVector, xold::AbstractVector): Computes the change in log-density between new and old coordinates.
- delta_logdensity(psi::YourWavefunctionType, xnew, xold::AbstractVector, dim::Int): Computes the change in log-density when only the coordinate at index dim is changed.
"""
abstract type AbstractWavefunction end

inputtype(::T) where {T <: AbstractWavefunction} = inputtype(T)
inputlength(::T) where {T <: AbstractWavefunction} = inputlength(T)
coordinate_projector(::T) where {T <: AbstractWavefunction} = coordinate_projector(T)
coordinate_projector(::Type{<:AbstractWavefunction}) = Base.identity
coordinate_transformer(::T) where {T <: AbstractWavefunction} = coordinate_transformer(T)
coordinate_transformer(::Type{<:AbstractWavefunction}) = Base.identity

@inline function _abs2(psi::AbstractWavefunction, x)
    return Base.abs2(psi(x))
end
@inline function logdensity(psi::AbstractWavefunction, x)
    return log(_abs2(psi, x))
end
@inline function delta_logdensity(psi::AbstractWavefunction, xnew::AbstractVector, xold::AbstractVector)
    return logdensity(psi, xnew) - logdensity(psi, xold)
end
@inline function delta_logdensity(psi::AbstractWavefunction, xnew, xold::AbstractVector, dim::Int)
    xnew_vec = copy(xold)
    xnew_vec[dim] = xnew
    return logdensity(psi, xnew_vec) - logdensity(psi, xold)
end

struct ProductWavefunction{D <: Tuple} <: AbstractWavefunction
    wavefunctions::D 
end
@inline function _abs2(psi::ProductWavefunction, x)
    result = _abs2(psi.wavefunctions[1], x)
    for wf in psi.wavefunctions[2:end]
        result *= _abs2(wf, x)
    end
    return result
end
@inline function logdensity(psi::ProductWavefunction, x)
    result = logdensity(psi.wavefunctions[1], x)
    for wf in psi.wavefunctions[2:end]
        result += logdensity(wf, x)
    end
    return result
end
@inline function delta_logdensity(psi::ProductWavefunction, xnew::AbstractVector, xold::AbstractVector)
    result = delta_logdensity(psi.wavefunctions[1], xnew, xold)
    for wf in psi.wavefunctions[2:end]
        result += delta_logdensity(wf, xnew, xold)
    end
    return result
end
@inline function delta_logdensity(psi::ProductWavefunction, xnew, xold::AbstractVector, dim::Int)
    result = delta_logdensity(psi.wavefunctions[1], xnew, xold, dim)
    for wf in psi.wavefunctions[2:end]
        result += delta_logdensity(wf, xnew, xold, dim)
    end
    return result
end
