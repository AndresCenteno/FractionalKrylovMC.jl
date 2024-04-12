module FractionalKrylovMC

# Write your package code here.
using LinearAlgebra
using AlphaStableDistributions, Distributions, Statistics
using Expokit, MittagLeffler
using HCubature
using ForwardDiff

include("types.jl")
include("utils.jl")
include("random_utils.jl")
include("solvers/MC.jl"); include("solvers/eigendecomposition.jl"); include("solvers/quadrature.jl")

export MittagLefflerProblem, FracExpProblem, MittagLefflerRand, FracExpRand
# abstract types to support different kinds of solutions
export MatVecSolver, MatVecSolution

# solving
export solve
export MCSolver
export MCSolverSaveSamples
export MCSolverExpokit
export QuadSolver
export EigenSolver

# utils 
export SpectralKernel, SpectralKernelRNG
export generate_times
export within_region
export create_random_problem
export relative_error

# random utils
export stblrnd
export stblrndsub

end
