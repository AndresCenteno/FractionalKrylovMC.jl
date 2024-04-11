
"""
Spectral kernel for the Mittag-Leffler function
"""
SpectralKernel(α,λ,r) = (λ*r^(α-1)*sin(α*π)) / (pi*(r^(2α)+2*λ*r^α*cos(α*π)+λ^2))

"""
RNG for γ-stable subordinator
"""
function γStableRNG(U::T,V::T,γ::T,t::T) where T<:Real
    α = γ
    """
    Computer Simulation of L6vy cz-Stable Variables and Processes
    Aleksander Weron * and Rafat Weron
    """
    S_α = (1 + tam)
end

"""
Generate ordered times for the quadrature of the Fractional Mittag-Leffler
and for the quadrature of the Fractional Exponential
"""
function generate_times(problem::MittagLefflerProblem{T},nsims::Int) where T<:Real
    λ, α, μ = problem.t^problem.α, problem.α, problem.μ
    LevySubordinator = AlphaStable(μ,1,cos(μ*pi/2)^(1/μ),0)
    SpectralKernelRNG(u) = abs(tan((α*pi)*u + atan(1/tan(pi*α)))*λ*sin(pi*α)-λ*cos(pi*α))^(1/α)
    times = SpectralKernelRNG.(rand(nsims)).^(1/μ).*rand(LevySubordinator,nsims)
    sort(filter(!isnan,times))
end

function generate_times(problem::FracExpProblem{T},nsims::Int;gradient=false) where T<:Real
    t, γ = problem.t, problem.γ
    if gradient == false
        LevySubordinator = AlphaStable(γ,1,cos(γ*pi/2)^(1/γ),0)
        times = t^(1/γ).*rand(LevySubordinator,nsims)
        return sort(filter(!isnan,times))
    end
    Ustream, Vstream = rand(nsims), rand(nsims)
    
end

"""
#TODO: why should the error follow the CLT? Well, if we take the mean, it actually should no?
Checks whether the MCSolution has the EigenSolution within its
(1-p) confidence region
https://en.wikipedia.org/wiki/Confidence_region
Using that if X∼N(μ,Σ)∈R^k then (X-μ)'*Σ^{-1}*(X-μ) follows a Χ^2 distribution
with k degrees of freedom.
See: https://stats.stackexchange.com/a/61993/404853
"""
function within_region(MCSol::MittagLefflerMCSolution,
    EigenSol::MittagLefflerEigenSolution,p::T) where T<:Real
    if MCSol.uT isa Vector
        error("Must provide matrix of samples, not just mean")
    end
    X = mean(MCSol.uT,dims=2)
    μ = EigenSol.uT
    Σ = cov(MCSol.uT,dims=2)
    k = length(μ); n = size(X,2) # degrees of freedom
    hotelling_sample = n*(X-μ)'*(Σ\(X-μ))
    confidence_interval = quantile(FDist(k,n),1-p)*k*(n-1)/(n-k)
    return hotelling_sample[1] < confidence_interval
end

struct MittagLefflerRand end
struct FracExpRand end

function create_random_problem(α::T,γ::T,k::Int,::MittagLefflerRand) where T<:Real
    A = rand(k,k); Q = Matrix(qr(A).Q); Λ = rand(k)
    t = rand()*2
    u0 = rand(k)
    problem = MittagLefflerProblem(Λ,Q,u0,t,α,γ)
    problem
end

function create_random_problem(γ::T,k::Int,::FracExpRand) where T<:Real
    A = rand(k,k); Q = Matrix(qr(A).Q); Λ = rand(k)
    t = rand()*2
    u0 = rand(k)
    problem = FracExpProblem(Λ,Q,u0,t,γ)
    problem
end


#TODO:
# Create struct called synthetic solution in which we have a positive semidefinite
# matrix created by ourselves and the orthogonal eigendecomposition is stored

relative_error(true_val,est_val) = norm(est_val-true_val)/norm(true_val)