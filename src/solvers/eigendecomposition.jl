"""
This is not really that useful hehe to be honest, we are realisticly never gonna
call this program
"""

struct EigenSolver <: MatVecSolver end

function solve(problem::MittagLefflerProblem{T}, ::EigenSolver
    )::MatVecSolution where T<:Real
    Λ, Χ = eigen(problem.A)
    α, γ = problem.α, problem.γ
    modified_spectrum = mittleff.(α,-Λ.^γ*problem.t^α)
    uT = Χ*diagm(modified_spectrum)*(Χ\problem.u0)
    return MittagLefflerEigenSolution(uT)
end