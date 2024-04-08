struct MCSolver <: MittagLefflerSolver end

function solve(problem::MittagLefflerProblem{T},nsims::Int,
    cutoff::T,::MCSolver)::MittagLefflerMCSolution where T<:Real
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

struct MCSolverSaveSamples <: MittagLefflerSolver end

function solve(problem::MittagLefflerProblem{T},nsims::Int,
    cutoff::T,::MCSolverSaveSamples)::MittagLefflerMCSolution where T<:Real
    uT_samples = zeros(problem.N,nsims)
    times = filter(x->x<=cutoff,generate_times(problem,nsims)) # times for quadrature
    jumps = [times[1];times[2:end]-times[1:end-1]]
    uT_samples[:,1] = exp(-problem.Anμ*jumps[1])*problem.u0
    for i in 2:length(jumps)
        uT_samples[:,i] = exp(-problem.Anμ*jumps[i])*uT_samples[:,i-1]
    end
    MittagLefflerMCSolution(uT_samples,times,nsims,cutoff)
end

struct MCSolverExpokit <: MittagLefflerSolver end

function solve(problem::MittagLefflerProblem{T},nsims::Int,
    cutoff::T,KrylovSize::Int,::MCSolverExpokit)::MittagLefflerMCSolution where T<:Real
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