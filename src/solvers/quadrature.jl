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
    uT = hcubature(u->exp_psd(TotalRNG([α,γ],u),-Anμ)*u0,[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    duTdp = (-Anμ)*hcubature(u->exp_psd(TotalRNG([α,γ],u),-Anμ)*u0*dTotalRNGdp([α;γ],u)',[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    return MittagLefflerSolution(uT,duTdp[:,1],duTdp[:,2])
end

struct QuadKrySolver <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T},KryDim::Int,::QuadKrySolver;atol=1e-5,rtol=1e-2) where T<:Real
    # unpack parameters
    Anμ = problem.Anμ; n = problem.nμ; γ = problem.γ; α = problem.α; t = problem.t; u0 = problem.u0
    # p is a vector of parameters p = [α;γ]
    TotalRNG(p,u) = SpectralKernelRNG(p[1],t^p[1],u[1])^(p[1]*n/p[2])*stblrndsub(p[2]/(p[1]*n),u[2],u[3])
    dTotalRNGdp(p,u) = ForwardDiff.gradient(p->TotalRNG(p,u),p)
    normu0 = norm(u0); u00 = u0/normu0;
    V, H = my_arnoldi(Anμ,KryDim,u00)
    e1 = zeros(size(H,1)); e1[1] = 1
    uT = hcubature(u->normu0*expmv_psd(TotalRNG([α;γ],u),V,-H,KryDim=KryDim),[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    uT = normu0*hcubature(u->(V*exp_cutoff(TotalRNG([α;γ],u),-H))[:,1],[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    duTdp = normu0*hcubature(u->(V*(-H)*exp_cutoff(TotalRNG([α;γ],u),-H))[:,1]*dTotalRNGdp([α;γ],u)',[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]

    # duTdp = (-Anμ)*hcubature(u->normu0*expmv_psd(TotalRNG([α;γ],u),V,-H,KryDim=KryDim)*dTotalRNGdp([α;γ],u)',[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    # duTdp = V*(-H)*V'*hcubature(u->normu0*expmv_psd(TotalRNG([α;γ],u),V,-H,KryDim=KryDim)*dTotalRNGdp([α;γ],u)',[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    # duTdp = hcubature(u->normu0*V*(-H)*exp(-TotalRNG([α;γ],u)*H)*e1*dTotalRNGdp([α;γ],u)',[0;0;0],[1;1;1],atol=atol,rtol=rtol)[1]
    # wtf am I doing here, think

    
    return MittagLefflerSolution(uT,duTdp[:,1],duTdp[:,2])
end

struct ShitSolver <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T},::ShitSolver;Nt=100) where T<:Real
    # unpack parameters
    Anμ = problem.Anμ; n = problem.nμ; γ = problem.γ; α = problem.α; t = problem.t; u0 = problem.u0
    TotalRNG(p,u) = SpectralKernelRNG(p[1],t^p[1],u[1])^(p[1]*n/p[2])*stblrndsub(p[2]/(p[1]*n),u[2],u[3])
    dTotalRNGdp(p,u) = ForwardDiff.gradient(p->TotalRNG(p,u),p)
    uT = shit_cuadrature_3D(u->exp_psd(TotalRNG([α,γ],u),-Anμ)*u0;Nt=Nt)
    duTdp = (-Anμ)*shit_cuadrature_3D(u->exp_psd(TotalRNG([α,γ],u),-Anμ)*u0*dTotalRNGdp([α;γ],u)')
    return MittagLefflerSolution(uT,duTdp[:,1],duTdp[:,2])
end