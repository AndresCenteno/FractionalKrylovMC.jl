"""
Here we are doing quadrature rules but instead of generating random numbers we are integrating
over [0,1]^3
"""

struct QuadSolver <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T},::QuadSolver;atol=1e-5,rtol=1e-2) where T<:Real
    # unpack parameters
    Anμ = problem.Anμ; n = problem.nμ; γ = problem.γ; α = problem.α; t = problem.t; u0 = problem.u0
    TotalRNG(p,u) = SpectralKernelRNG(p[1],t^p[1],u[1])^(p[1]*n/p[2])*stblrndsub(p[2]/(p[1]*n),u[2],u[3])
    dTotalRNGdp(p,u) = ForwardDiff.gradient(p->TotalRNG(p,u),p)
    uT = hcubature(u->exp(-Anμ*TotalRNG([α,γ],u))*u0,[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    duTdp = (-Anμ)*hcubature(u->exp(-Anμ*TotalRNG([α;γ],u))*u0*dTotalRNGdp([α;γ],u)',[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    return MittagLefflerSolution(uT,duTdp[:,1],duTdp[:,2])
end

struct QuadKrySolver <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T},KryDim::Int,::QuadKrySolver;atol=1e-5,rtol=1e-2) where T<:Real
    # unpack parameters
    Anμ = problem.Anμ; n = problem.nμ; γ = problem.γ; α = problem.α; t = problem.t; u0 = problem.u0
    TotalRNG(p,u) = SpectralKernelRNG(p[1],t^p[1],u[1])^(p[1]*n/p[2])*stblrndsub(p[2]/(p[1]*n),u[2],u[3])
    dTotalRNGdp(p,u) = ForwardDiff.gradient(p->TotalRNG(p,u),p)
    uT = hcubature(u->expmv(TotalRNG([α;γ],u),-Anμ,u0,m=KryDim),[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    duTdp = (-Anμ)*hcubature(u->expmv(TotalRNG([α;γ],u),-Anμ,u0,m=KryDim)*dTotalRNGdp([α;γ],u)',[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    return MittagLefflerSolution(uT,duTdp[:,1],duTdp[:,2])
end