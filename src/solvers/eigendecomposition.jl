struct EigenSolver <: MittagLefflerSolver end

function solve(problem::MittagLefflerProblem{T}, ::EigenSolver
    )::MittagLefflerEigenSolution where T<:Real
    Λ, Χ = eigen(problem.A)
    α, γ = problem.α, problem.γ
    modified_spectrum = mittleff.(α,-Λ.^γ*problem.t^α)
    uT = Χ*diagm(modified_spectrum)*(Χ\problem.u0)
    return MittagLefflerEigenSolution(uT)
end