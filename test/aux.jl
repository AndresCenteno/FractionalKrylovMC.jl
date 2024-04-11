using FractionalKrylovMC
using Test
# do not want to export this from FractionalKrylovMC within that namespace
using LinearAlgebra: cholesky
using Distributions: FDist
using AlphaStableDistributions
using Statistics

    test_passed = 0
    experiments = Int(2)
    nsims = Int(1e5)
    for _ in 1:experiments
        λ = rand()*2
        t = rand()
        γ = 0.1 + rand()*0.9
        aux = λ*t^(1/γ)
        τ = AlphaStable(2*γ,1,cos(γ*pi/2)^(1/2*γ),0)
        samples = zeros(nsims)
        for sim in 1:nsims
            samples[sim] = exp(-aux*rand(τ))
        end
        test_passed += abs(exp(-λ*t^γ)-mean(samples)) <= std(samples)*1.96/sqrt(nsims)
    end
    @show test_passed, experiments

    exp(-λ*t^γ)
    std(samples)*1.96/sqrt(nsims)
    mean(samples)
    τ = AlphaStable(γ,1,cos(γ*pi/2)^(1/γ),0)

    rand(τ)