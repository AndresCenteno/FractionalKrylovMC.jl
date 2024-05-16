using FractionalKrylovMC

problem = create_random_problem(rand()*0.9 + 0.1,rand()*0.9 + 0.1,10,MittagLefflerRand())
uT, list_of_times = solve(problem,IntegratorSolver(); Nt = 200)
problem.solution.uT
uT_quad = solve(problem,QuadSolver())
uT_quad.uT
problem.u0
println(list_of_times)