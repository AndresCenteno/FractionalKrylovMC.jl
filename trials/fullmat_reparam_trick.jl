using FractionalKrylovMC

problem = create_random_problem(0.3,0.5,50,MittagLefflerRand())

quad_solution = solve(problem,QuadSolver())

quadkry_solution = solve(problem,4,QuadKrySolver())

relative_error(problem.solution,quad_solution)
relative_error(problem.solution,quadkry_solution)
quad_solution.duTdα
problem.solution.duTdα