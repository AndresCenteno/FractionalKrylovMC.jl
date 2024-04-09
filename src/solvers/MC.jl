struct MCSolver <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T},nsims::Int,
    cutoff::T,::MCSolver)::MatVecSolution where T<:Real
    uT = zeros(problem.N)
    times = filter(x->x<=cutoff,generate_times(problem,nsims)) # times for quadrature
    jumps = [times[1];times[2:end]-times[1:end-1]]
    v = copy(problem.u0)
    for i in eachindex(jumps)
        v = exp(-problem.Anμ*jumps[i])*v
        uT += v/nsims # times beyond cutoff are considered steady state 0
    end
    MittagLefflerMCSolution(uT,times,nsims,cutoff)
end

function solve(problem::FracExpProblem{T},nsims::Int,
    cutoff::T,::MCSolver)::MatVecSolution where T<:Real
    uT = zeros(problem.N)
    times = filter(x->x<=cutoff,generate_times(problem,nsims)) # times for quadrature
    jumps = [times[1];times[2:end]-times[1:end-1]]
    v = copy(problem.u0)
    for i in eachindex(jumps)
        v = exp(-problem.A*jumps[i])*v
        uT += v/nsims # times beyond cutoff are considered steady state 0
    end
    FracExpMCSolution(uT,times,nsims,cutoff)
end

struct MCSolverSaveSamples <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T},nsims::Int,
    cutoff::T,::MCSolverSaveSamples)::MatVecSolution where T<:Real
    uT_samples = zeros(problem.N,nsims)
    times = filter(x->x<=cutoff,generate_times(problem,nsims)) # times for quadrature
    jumps = [times[1];times[2:end]-times[1:end-1]]
    uT_samples[:,1] = exp(-problem.Anμ*jumps[1])*problem.u0
    for i in 2:length(jumps)
        uT_samples[:,i] = exp(-problem.Anμ*jumps[i])*uT_samples[:,i-1]
    end
    MittagLefflerMCSolution(uT_samples,times,nsims,cutoff)
end

function solve(problem::FracExpProblem{T},nsims::Int,
    cutoff::T,::MCSolverSaveSamples)::MatVecSolution where T<:Real
    uT_samples = zeros(problem.N,nsims)
    times = filter(x->x<=cutoff,generate_times(problem,nsims)) # times for quadrature
    jumps = [times[1];times[2:end]-times[1:end-1]]
    uT_samples[:,1] = exp(-problem.Anμ*jumps[1])*problem.u0
    for i in 2:length(jumps)
        uT_samples[:,i] = exp(-problem.A*jumps[i])*uT_samples[:,i-1]
    end
    FracExpMCSolution(uT_samples,times,nsims,cutoff)
end

struct MCSolverExpokit <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T},nsims::Int,
    cutoff::T,KrylovSize::Int,::MCSolverExpokit)::MatVecSolution where T<:Real
    uT = zeros(problem.N)
    times = filter(x->x<=cutoff,generate_times(problem,nsims)) # times for quadrature
    jumps = [times[1];times[2:end]-times[1:end-1]]
    v = copy(problem.u0)
    for i in eachindex(jumps)
        v = expmv(jumps[i],-problem.Anμ,v,m=KrylovSize)
        uT += v/nsims # times beyond cutoff are considered steady state 0
    end
    MittagLefflerMCSolution(uT,times,nsims,cutoff)
end

function solve(problem::FracExpProblem{T},nsims::Int,
    cutoff::T,KrylovSize::Int,::MCSolverExpokit)::MatVecSolution where T<:Real
    uT = zeros(problem.N)
    times = filter(x->x<=cutoff,generate_times(problem,nsims)) # times for quadrature
    jumps = [times[1];times[2:end]-times[1:end-1]]
    v = copy(problem.u0)
    for i in eachindex(jumps)
        v = expmv(jumps[i],-problem.A,v,m=KrylovSize)
        uT += v/nsims # times beyond cutoff are considered steady state 0
    end
    FracExpMCSolution(uT,times,nsims,cutoff)
end