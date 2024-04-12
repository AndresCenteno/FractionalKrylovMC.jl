using FractionalKrylovMC

problem = create_random_problem(0.3,0.5,50,MittagLefflerRand())

quad_solution = solve(problem,QuadSolver(),atol=1e-7,rtol=1e-4)

quadkry_solution = solve(problem,4,QuadKrySolver(),atol=1e-7,rtol=1e-4)

relative_error(problem.solution,quad_solution)
# 6.393421459972954e-6
# 0.035779317152599306
# 0.03216800727458906
relative_error(problem.solution,quadkry_solution)
# 7.640089340742012e-6
# 0.04963438749727011
# 0.6539085724845574