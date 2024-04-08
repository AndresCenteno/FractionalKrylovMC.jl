module FractionalKrylovMC

# Write your package code here.
using LinearAlgebra
using AlphaStableDistributions, Distributions, Statistics
using ExponentialUtilities, MittagLeffler

include("types.jl")
include("utils.jl")
include("solvers/MC.jl"); include("solvers/eigendecomposition.jl")

export MittagLefflerProblem
# abstract types to support different kinds of solutions
export MittagLefflerSolver
export MittagLefflerSolution
export MittagLefflerMCSolution

# solving
export solve
export MCSolver
export MCSolverSaveSamples
export EigenSolver

# utils
export SpectralKernel
export generate_times
export within_region
export create_random_problem

end
