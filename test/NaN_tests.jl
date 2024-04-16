# checking where the NaNs come from

using FractionalKrylovMC

flag = 1
k = 30
while true
    @show flag
    flag += 1
    global problem = create_random_problem(rand()*0.9+0.1,rand()*0.9+0.1,k,MittagLefflerRand())
    try
        sol = solve(problem,10,QuadKrySolver())
    catch
        break
    end
end

@run sol = solve(problem,10,QuadKrySolver())

