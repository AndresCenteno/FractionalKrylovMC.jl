using FractionalKrylovMC

problem = create_random_problem(rand()*0.9 + 0.1,rand()*0.9 + 0.1,10,MittagLefflerRand())
qsolution = solve(problem,ShitSolver();Nt=50)
qsolution.uT
problem.solution.uT
# https://www.value-at-risk.net/numerical-integration-multiple-dimensions/