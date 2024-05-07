using FractionalKrylovMC

problem = create_random_problem(rand()*0.9 + 0.1,rand()*0.9 + 0.1,10,MittagLefflerRand())

problem.solution.uT
α -> E_α(-A^γ*t^α)*u_0
problem.solution.duTdα
problem.solution.duTdγ
problem.A
