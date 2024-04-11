# I'm doing on .jl instead of .ipynb because no one in the Julia community uses .ipynb don't know why

using FractionalKrylovMC

# we define a problem by
γ = 0.7 # space exponent
k = 200 # size of matrix
problem = create_random_problem(γ,k,FracExpRand())

# the theoretical solution is already in the problem 
problem.uT_spectral

# we invoke the MCSolver
nsims = Int(1e4)
cutoff = 1e6
MCSol = solve(problem,nsims,cutoff,MCSolver())

relative_error(problem.uT_spectral,MCSol.uT)

# and then the Krylov MCSolver
KrySubspaceDim = 5
MCKrySol = solve(problem,nsims,cutoff,KrySubspaceDim,MCSolverExpokit())

relative_error(problem.uT_spectral,MCKrySol.uT)

#TODO: convergence of error