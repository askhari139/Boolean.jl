module Boolean

using Pkg
using CSV
using DataFrames
using SparseArrays
using LinearAlgebra
using Pipe
using ProgressMeter
using Base.Threads
using Random
using Combinatorics
using JSON
using PyCall
import IterTools: product
# using Plots

# include("dependencies.jl")
include("utils.jl")
include("async_update.jl")
include("multiLevel_shubham.jl")
include("CSB.jl")
include("async_non_matrix.jl")
include("customFunctions.jl")
include("oddLevel.jl")
include("logicalRules.jl")
include("logicalSim.jl")
include("bmodel.jl")
include("boolToTopo.jl")
include("coherence.jl")
include("stochasticSim.jl")


# Determine the absolute path to the Python script's directory
script_dir = joinpath(@__DIR__)  # @__DIR__ gives path to the current file (src/)

# Add that directory to Python's import path
push!(PyVector(pyimport("sys")."path"), "/Users/kishorehari/Desktop/Boolean.jl/src")

# Import your script (assumes it's named `script.py`)
logicRules = pyimport("logicRules")


export bmodel_reps,
    edgeWeightPert,
    weightedTopoSim,
    contWeightPert,
    getSSListRand,
    scanNodeTurnOff,
    getNodes,
    coherence,
    coherenceAllNode

end # module Boolean
