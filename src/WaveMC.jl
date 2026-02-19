"""
Allocation-free Quantum Monte Carlo for arbitrary wave functions 
and observables in Julia, built on the Carlo.jl framework.
"""
module WaveMC


export ComplexNormal
export AbstractWavefunctionMC, WavefunctionMC
export BufferedObservables
export Histogram, plot_histogram

using Carlo
using Printf

using Distributions
using Random
using HDF5

using ThreadsX
using ChunkSplitters
using CairoMakie
using DataFrames

include("wavefunction/wavefunction.jl")
include("acceptance_rate/acceptance_rate_adapter.jl")
include("buffer/buffer.jl")
include("coordinate_update/coordinate_update.jl")
include("state.jl")
include("AbstractWavefunctionMC.jl")
include("acceptance_rate/acceptance_rate.jl")
include("observables/abstractobservables.jl")
include("complex_normal_distribution.jl")
include("WavefunctionMC.jl")

include("postprocess/binning.jl")
include("postprocess/histogram.jl")
include("postprocess/plot.jl")

end
