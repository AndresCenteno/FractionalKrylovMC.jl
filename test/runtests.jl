using FractionalKrylovMC
using Test
# do not want to export this from FractionalKrylovMC within that namespace
using LinearAlgebra: cholesky
using Distributions: FDist
using Statistics

# this test does not pass because of covariancee matrix not being invertible
@testset "FractionalKrylovMC.jl" begin
    p = 0.1
    k = 20
    experiments = Int(25)
    nsims = Int(1e6)
    cutoff = 1e7
    test_passed = 0
    for _ in 1:experiments
        # generate random parameters
        problem = create_random_problem(k)
        MCSol = solve(problem,nsims,cutoff,MCSolverSaveSamples())
        EigenSol = solve(problem,EigenSolver())
        test_passed += within_region(MCSol,EigenSol,p)
    end
    @show test_passed, experiments
    @test test_passed >= (1-p)*experiments
end

@testset "ConfidenceIntervals" begin
    """
    Test for 1e2 experiments drawing 1e6 samples from N(μ,Σ) with Σ unknown
    if we pass the tests for
    """
    p = 0.05
    k = 30
    experiments = Int(3e2)
    samples_per_experiment = Int(1e5); n = samples_per_experiment
    test_passed = 0
    confidence_interval = quantile(FDist(k,nsims),1-p)*k*(nsims-1)/(nsims-k)
    for _ in 1:experiments
        # mean μ does not matter at all here
        Σ = rand(k,k); Σ = Σ*Σ'
        σ = cholesky(Σ).L
        samples = σ*randn(k,n)
        X = mean(samples,dims=2)
        Σsampled = cov(samples,dims=2) # = (samples .- X)*(samples .- X)'/(n-1)
        hotelling_sample = n*X'*(Σsampled\X)
        test_passed += hotelling_sample[1] <= confidence_interval
    end
    @show test_passed, experiments
    @test test_passed >= 0.95*experiments
end 

# @testset "IntegralFormulas" begin
#     """
#     I pretend to test the case γ > α
#     """
#     γ = 0.8; α = 0.3
    
# end