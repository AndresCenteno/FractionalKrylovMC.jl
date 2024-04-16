using FractionalKrylovMC
using Test
# do not want to export this from FractionalKrylovMC within that namespace
using LinearAlgebra: cholesky, qr, diagm
# using Distributions: FDist
using Statistics

@testset "No Krylov proyection" begin
    k = 20 # size of matrix
    experiments = Int(5)
    for _ in 1:experiments
        # generate random parameters
        problem = create_random_problem(rand()*0.9+0.1,rand()*0.9+0.1,k,MittagLefflerRand())
        MCSol = solve(problem,QuadSolver(),atol=1e-8,rtol=1e-5)
        @test all(relative_error(problem.solution,MCSol).<=1e-2)
    end
    # one test fails :P, it's not really a lot, 5 experiments is almost 3 min with this much rtol
    # but keep in mind that
end

@testset "Moderate k=100 -> k=20 Krylov proyection" begin
    k = 100 # size of matrix
    kkry = 20 # Krylov downsize
    experiments = Int(3)
    for _ in 1:experiments
        # generate random parameters
        problem = create_random_problem(rand()*0.9+0.1,rand()*0.9+0.1,k,MittagLefflerRand())
        MCSol = solve(problem,kkry,QuadKrySolver(),atol=1e-8,rtol=1e-5)
        @test all(relative_error(problem.solution,MCSol).<=1e-2)
    end
end

# no puede ser que la formula integral no me funcione madre mia
@testset "IntegralFormulas" begin
    """
    Checking in 1D whether exp(-t*λ^γ) = exp(-(t^(1/γ)*λ)^γ)=∫_0^∞ exp(-t^(1/γ)*λ*τ)p(γ,τ)dτ
    using the CLT
    """
    test_passed = 0;
    tests = 200;
    nsims = Int(1e4);
    error_vec = zeros(tests,1);
    for test = 1:tests
        alpha = rand()*0.9 + 0.1;
        skewness = cos(alpha*pi/2)^(1/alpha);
        samples = zeros(nsims,1);
        lambda = rand(); t = rand();
        lambdatalpha = t^(1/alpha)*lambda;
        for sim = 1:nsims
            tau = stblsubrnd(alpha);
            samples[sim] = exp(-tau*lambdatalpha);
        end
        error_vec[test] = abs(mean(samples)-exp.(-t*lambda^alpha))/abs(exp.(-t*lambda^alpha));
        boolean = abs(mean(samples)-exp.(-t*lambda^alpha))<= 1.96*std(samples)/sqrt(nsims);
        test_passed = test_passed + boolean;
    end
    @test test_passed >= 0.94*tests
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