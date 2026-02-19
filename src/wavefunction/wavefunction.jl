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
    return prod(_abs2(wf, x) for wf in psi.wavefunctions)
end
@inline function logdensity(psi::ProductWavefunction, x)
    return sum(logdensity(wf, x) for wf in psi.wavefunctions)
end
@inline function delta_logdensity(psi::ProductWavefunction, xnew::AbstractVector, xold::AbstractVector)
    return sum(delta_logdensity(wf, xnew, xold) for wf in psi.wavefunctions)
end
@inline function delta_logdensity(psi::ProductWavefunction, xnew, xold::AbstractVector, dim::Int)
    return sum(delta_logdensity(wf, xnew, xold, dim) for wf in psi.wavefunctions)
end

struct SemiDirectProductWavefunction{D <: Tuple, T <: Tuple, N, S, coordinate_proj, coordinate_trans} <: AbstractWavefunction
    wavefunctions::D
    coordinates::T ## Each wavefunction takes the positions of its coordinates as input, e.g. coordinates = (1:N-1, N:N)
    coordinate_dict::Dict{Int, Tuple{Int,Int}} ## Maps coordinate index to the wavefunction it belongs to and the dim inside that subarray, e.g. Dict(1=>1, 2=>1, 3=>2) for coordinates = (1:2, 3:3)
    function SemiDirectProductWavefunction(wavefunctions::D, coordinates::T) where {D <: Tuple, T <: Tuple}
        coordinate_dict = Dict{Int, Tuple{Int,Int}}()
        for (i, coords) in enumerate(coordinates)
            for (j,c) in enumerate(coords)
                coordinate_dict[c] = (i,j)
            end
        end
        N = sum(length.(coordinates))
        S = only(unique(inputtype.(wavefunctions)))
        coordinate_proj = coordinate_projector(wavefunctions[1])
        coordinate_trans= coordinate_transformer(wavefunctions[1])
        return new{D, T, N, S, coordinate_proj, coordinate_trans}(wavefunctions, coordinates, coordinate_dict)
    end
end
inputtype(::SemiDirectProductWavefunction{D, T, N, S}) where {D, T, N, S} = S
inputlength(::SemiDirectProductWavefunction{D, T, N, S}) where {D, T, N, S} = N
coordinate_projector(psi::SemiDirectProductWavefunction{D, T, N, S, coordinate_proj, coordinate_trans}) where {D, T, N, S, coordinate_proj, coordinate_trans} = coordinate_proj
coordinate_transformer(psi::SemiDirectProductWavefunction{D, T, N, S, coordinate_proj, coordinate_trans}) where {D, T, N, S, coordinate_proj, coordinate_trans} = coordinate_trans



@inline function _abs2(psi::SemiDirectProductWavefunction, x)
    return prod(_abs2(wf, view(x, psi.coordinates[i])) for (i,wf) in enumerate(psi.wavefunctions))
end
@inline function logdensity(psi::SemiDirectProductWavefunction, x)
    return sum(logdensity(wf, view(x, psi.coordinates[i])) for (i,wf) in enumerate(psi.wavefunctions))
end
@inline function delta_logdensity(psi::SemiDirectProductWavefunction, xnew::AbstractVector, xold::AbstractVector)
    return sum(delta_logdensity(wf, view(xnew, psi.coordinates[i]), view(xold, psi.coordinates[i])) for (i,wf) in enumerate(psi.wavefunctions))
end
@inline function delta_logdensity(psi::SemiDirectProductWavefunction, xnew, xold::AbstractVector, dim::Int)
    ## First, we need to find out in which wavefunction the coordinate is changed! All the others do not contribute here!
    ind, dim_subarray = psi.coordinate_dict[dim]
    return delta_logdensity(psi.wavefunctions[ind], xnew, view(xold, psi.coordinates[ind]), dim_subarray)
end
