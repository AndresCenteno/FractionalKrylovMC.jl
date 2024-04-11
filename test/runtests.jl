using FractionalKrylovMC
using Test
# do not want to export this from FractionalKrylovMC within that namespace
using LinearAlgebra: cholesky, qr, diagm
# using Distributions: FDist
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
    λtγ = λ*t^(1/γ)
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

# no puede ser que la formula integral no me funcione madre mia
@testset "IntegralFormulas" begin
    """
    Checking in 1D whether exp(-t*λ^γ) = exp(-(t^(1/γ)*λ)^γ)=∫_0^∞ exp(-t^(1/γ)*λ*τ)p(γ,τ)dτ
    using the CLT
    """
    # test_passed = 0;
    # tests = 300;
    # nsims = Int(1e5);
    # error_vec = zeros(tests,1);
    # for test = 1:tests
    #     alpha = rand()*0.9 + 0.1;
    #     skewness = cos(alpha*pi/2)^(1/alpha);
    #     samples = zeros(nsims,1);
    #     lambda = rand(); t = rand();
    #     lambdatalpha = t^(1/alpha)*lambda;
    #     for sim = 1:nsims
    #         tau = stblsubrnd(alpha);
    #         samples[sim] = exp(-tau*lambdatalpha);
    #     end
    #     error_vec[test] = abs(mean(samples)-exp.(-t*lambda^alpha))/abs(exp.(-t*lambda^alpha));
    #     boolean = abs(mean(samples)-exp.(-t*lambda^alpha))<= 1.96*std(samples)/sqrt(nsims);
    #     test_passed = test_passed + boolean;
    # end
    # @test test_passed >= 0.94*tests
    """
    Testing same as above but for one matrix
    """
    k = 20; γ = 0.1 + rand()*0.9
    λs = rand(k); A = rand(k,k); xs = Matrix(qr(A).Q); u0 = rand(k);
    reconstructedM = xs*diagm(λs)*xs'
    uT_spectral = xs*diagm(-λs.^γ)*xs'*u0
    uT_sims = zeros(k); nsims = Int(1e6)
    random_times = stblrndsub(γ,nsims)
    sort(filter())
    for sim=1:nsims
        uT_sims += exp(-reconstructedM*random_times[sim])*u0/nsims
    end
    @test isapprox(uT_sims, uT_spectral; rtol = 1e-2)
end

@testset "Type stability" begin
    """
    Testing small stuff about the stblsubrnd
    """
    alpha = BigFloat(rand()*0.9 + 0.1)
    @test stblsubrnd(alpha) isa BigFloat
    @test stblsubrnd(alpha,2) isa Vector{BigFloat}
    @test stblsubrnd(Float64(alpha),rand(Float64,3),rand(Float64,3)) isa Vector{Float64}
end