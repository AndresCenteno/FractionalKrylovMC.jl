module FractionalKrylovMC

# Write your package code here.
using LinearAlgebra
using LinearAlgebra: norm2
using AlphaStableDistributions, Distributions, Statistics
using Expokit, MittagLeffler
using HCubature
using ForwardDiff

include("types.jl")
include("utils.jl")
include("random_utils.jl")
include("solvers/MC.jl"); include("solvers/eigendecomposition.jl"); include("solvers/quadrature.jl")

export MittagLefflerProblem, FracExpProblem, MittagLefflerRand, FracExpRand, MittagLefflerSolution
# abstract types to support different kinds of solutions
export MatVecSolver, MatVecSolution

# solving
export solve
export MCSolver
export MCSolverSaveSamples
export MCSolverExpokit
export QuadSolver, QuadKrySolver, ShitSolver
export EigenSolver

# utils 
export SpectralKernel, SpectralKernelRNG
export generate_times
# export within_region # it doesn't make any sense to explore convergence of CLT when we can just
# integrate deterministically over the RNG stream
export create_random_problem
export relative_error

# random utils
export stblrnd
export stblrndsub
export SpectralKernelRNG, SpectralKernel

end
