using FractionalKrylovMC
using Test

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
    experiments = Int(1e2)
    samples_per_experiment = Int(1e6)
    test_passed = 0
    for _ in 1:experiments
        Σ = rand()
    end
end 