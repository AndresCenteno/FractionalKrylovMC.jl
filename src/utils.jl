
"""
Spectral kernel for the Mittag-Leffler function
"""
SpectralKernel(α,λ,r) = (λ*r^(α-1)*sin(α*π)) / (pi*(r^(2α)+2*λ*r^α*cos(α*π)+λ^2))
SpectralKernelRNG(α,λ,u) = abs(tan((α*pi)*u + atan(1/tan(pi*α)))*λ*sin(pi*α)-λ*cos(pi*α))^(1/α)

"""
Generate ordered times for the quadrature of the Fractional Mittag-Leffler
and for the quadrature of the Fractional Exponential
"""
function generate_times(problem::MittagLefflerProblem{T},nsims::Int) where T<:Real
    λ, α, μ = problem.t^problem.α, problem.α, problem.μ
    SpectralKernelRNG(u) = abs(tan((α*pi)*u + atan(1/tan(pi*α)))*λ*sin(pi*α)-λ*cos(pi*α))^(1/α)
    times = SpectralKernelRNG.(rand(nsims)).^(1/μ).*stblrndsub(μ,nsims)
    sort(filter(!isnan,times))
end

function generate_times(problem::FracExpProblem{T},nsims::Int;gradient=false) where T<:Real
    t, γ = problem.t, problem.γ
    if gradient == false
        times = t^(1/γ).*stblrndsub(μ,nsims)
        return sort(filter(!isnan,times))
    end    
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
    # A = Q*Λ*inv(Q)
    # f(A) = Q*f(Λ)*inv(Q)
    t = rand()*2
    u0 = rand(k)
    problem = FracExpProblem(Λ,Q,u0,t,γ)
    problem
end


#TODO:
# Create struct called synthetic solution in which we have a positive semidefinite
# matrix created by ourselves and the orthogonal eigendecomposition is stored

relative_error(true_val::T,est_val::T) where T<:AbstractVector = norm2(est_val-true_val)/norm2(true_val)

function relative_error(true_solution::MittagLefflerSolution{T},quad_solution::MittagLefflerSolution{T}) where T<:Real
    sol_err = zeros(3)
    sol_err[1] = relative_error(true_solution.uT,quad_solution.uT)
    sol_err[2] = relative_error(true_solution.duTdα,quad_solution.duTdα)
    sol_err[3] = relative_error(true_solution.duTdγ,quad_solution.duTdγ)
    sol_err
end

#### AVOIDING LAPACK ERRORS
function exp_psd(t::T,A::Matrix{T};cutoff::T=1e7) where T<:Real
    if t < cutoff
        return exp(A*t)
    end
    return fill!(similar(A),zero(T))
end

function expmv_psd(t::T,V::Matrix{T},H::Matrix{T};KryDim::Int,cutoff::T=1e7) where T<:Real
    if t < cutoff
        return (V*exp(H*t))[:,1]
    end
    return zeros(size(V,1))
end

function my_arnoldi(A::Matrix{T},KryDim::Int,v::Vector{T}) where T<:Real
    n = size(A,1);
    V = zeros(n,KryDim)
    H = zeros(KryDim,KryDim)
    
    V[:,1] = v
    for k=2:KryDim+1
        w = A*V[:,k-1]
        for j=1:k-1
            b = V[:,j]
            dotProd = (w'*b)
            w = w .- dotProd*b
            H[j,k-1] = dotProd
        end
        if k != KryDim + 1
            H[k,k-1] = norm(w)
            V[:,k] = w./norm(w)
        end
    end
    return V, H
end