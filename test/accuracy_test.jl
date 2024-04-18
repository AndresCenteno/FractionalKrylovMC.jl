using FractionalKrylovMC
using LinearAlgebra: norm2
using Statistics: mean, std
using DelimitedFiles

###################################
######## EXPERIMENT 1 #############
###################################

# fixing histograms, 3 histograms, constant size, random alpha and gamma, no Krylov
nsims = 10
samples = zeros(3,nsims)
k = 50
my_rtol = 1e-4; my_atol = 1e-4
# checking relative error
for sim = 1:nsims
    @show sim
    problem = create_random_problem(rand()*0.9 + 0.1,rand()*0.9 + 0.1,k,MittagLefflerRand())
    qsolution = solve(problem,QuadSolver(),atol=my_atol,rtol=my_rtol)
    samples[:,sim] = relative_error(problem.solution,qsolution)
end

for i=1:3
    println("$(mean((filter(!isnan,samples[i,:])))) ± $(std((filter(!isnan,samples[i,:]))))")
end

# with rtol = 1e-4, atol = 1e-4
# 6.4398226531423046e-6 ± 6.1938064777916645e-6
# 0.000487037793684316 ± 0.0015983058975240737
# 0.0002551078373981897 ± 0.0008680105437434714

length(filter(isnan,samples[1,:])) # 6/30 NaNs I don't know why

writedlm("test/accuracy_experiments/experiment1.csv",samples,',')

###################################
######## EXPERIMENT 2 #############
###################################

# same as experiment 1 but with Krylov

nsims = 3
k = 200
KryDimVec = 5
# hacer df

samples = zeros(3,nsims)
exponents = zeros(2,nsims)

my_rtol = 1e-5
# checking relative error
for sim = 1:nsims
    @show sim
    α = rand()*0.9 + 0.1; γ = rand()*0.9 + 0.1
    exponents[:,sim] = [α,γ]
    @show [α,γ]
    problem = create_random_problem(α,γ,k,MittagLefflerRand())
    qsolution = solve(problem,KryDimVec,QuadKrySolver(),rtol=my_rtol)
    samples[:,sim] = relative_error(problem.solution,qsolution)
end

for i=1:3
    println("$(mean((filter(!isnan,samples[i,:])))) ± $(std((filter(!isnan,samples[i,:]))))")
end

# NaN thing is solved
# 0.0050924031332660645 ± 0.007432084442461333
# 0.00996135216351648 ± 0.014169568713338587
# 0.016644908556880805 ± 0.03126378903737406
writedlm("test/accuracy_experiments/experiment3.csv",samples,',')

samples[:,:,1]