using Plots, FractionalKrylovMC
using LinearAlgebra: norm2
using DelimitedFiles
using Statistics: mean, std
using DataFrames

###################################
######## EXPERIMENT 1 #############
###################################

# fixing histograms, 3 histograms, constant size, random alpha and gamma, no Krylov
nsims = 30
samples = zeros(3,30)
k = 20
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

nsims = 4
k = 70
KryDimVec = 10:5:25
# hacer df

samples = zeros(3,nsims,length(KryDimVec))
exponents = zeros(2,nsims)

my_rtol = 1e-2
# checking relative error
for sim = 1:nsims
    @show sim
    α = rand()*0.3 + 0.7; γ = rand()*0.3 + 0.7
    exponents[:,sim] = [α,γ]
    @show [α,γ]
    problem = create_random_problem(α,γ,k,MittagLefflerRand())
    for i in eachindex(KryDimVec)
        @show i
        qsolution = solve(problem,KryDimVec[i],QuadKrySolver(),rtol=my_rtol)
        samples[:,sim,i] = relative_error(problem.solution,qsolution)
    end
end

# [α, γ] = [0.3302356643852082, 0.5253787346737683]
(α,γ) = [0.3302356643852082, 0.5253787346737683]

# I get NaNs with those numbers, my guess is: I get them at the rng, trying to see where NaNs are coming from
using Revise
using FractionalKrylovMC
SpectralKernelRNG.(α,1,range(0,1,10000)) # no NaNs here
n = ceil(γ/α)
μ = γ/(α*n)
Δt = 1e-2
times = Δt:Δt:(1-Δt)
rngres = [stblrndsub(μ,u,v) for u in times, v in times]
plot(times,times,rngres,st=:surface,camera=(40,60))
using ForwardDiff
drngresdα = [ForwardDiff.derivative(α->stblrndsub(γ/(α*2),u,v),α) for u in times, v in times]
drngresdγ = [ForwardDiff.derivative(γ->stblrndsub(γ/(α*2),u,v),γ) for u in times, v in times]

plot(times,times,drngresdα,st=:surface,camera=(40,60))
plot(times,times,drngresdγ,st=:surface,camera=(40,60))
totaldrngresdα = [ForwardDiff.derivative(α->SpectralKernel(α,2^α,w)*stblrndsub(γ/(α*2),u,v),α) for u in times, v in times, w in times]

filter(isnan,totaldrngresdα)
# maybe try to emulate by generating a bunch of problems and catching the nans

flag = 1
while flag > 0
    @show flag
    flag += 1
    global problem = create_random_problem(rand()*0.9+0.1,rand()*0.9+0.1,20,MittagLefflerRand())
    try
        solve(problem,QuadSolver(),rtol=1e-5)
    catch
        flag = 0
    end
end

using FractionalKrylovMC
flag = 1
while true
    @show flag
    flag += 1
    global problem = create_random_problem(rand()*0.9+0.1,rand()*0.9+0.1,20,MittagLefflerRand())
    sol = solve(problem,QuadSolver(),rtol=1e-5)
    if any(isnan.(sol.uT))
        break
    end
end
problem
problem.α
problem.γ
problem.N
solve(problem,QuadSolver())
solve(problem,20,QuadKrySolver())