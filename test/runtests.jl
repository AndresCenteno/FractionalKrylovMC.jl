using FractionalKrylovMC
using Test
# do not want to export this from FractionalKrylovMC within that namespace
using LinearAlgebra: cholesky
using Distributions: FDist

@testset "FractionalKrylovMC.jl" begin
    # Write your tests here.
    N = 20
    A = rand(N,N); A = A*A';
    α = 0.9; γ = 0.8
    t = 2.
    u0 = randn(N)
    problem = MittagLefflerProblem(A,u0,t,α,γ)
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
    confidence_interval = quantile(FDist(k,n),1-p)*k*(n-1)/(n-k)
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