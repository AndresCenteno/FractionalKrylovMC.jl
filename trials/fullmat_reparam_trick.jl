using FractionalKrylovMC

problem = create_random_problem(0.3,0.5,20,MittagLefflerRand())

quad_solution = solve(problem,QuadSolver())

relative_error(problem.solution,quad_solution)

quad_solution.duTdα
problem.solution.duTdα