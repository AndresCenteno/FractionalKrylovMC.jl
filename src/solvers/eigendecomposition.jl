struct EigenSolver <: MittagLefflerSolver end

function solve(problem::MittagLefflerProblem{T}, ::EigenSolver
    )::MittagLefflerEigenSolution where T<:Real
    Λ, Χ = eigen(problem.A)
    modified_spectrum = mittleff.(problem.α,-Λ.^problem.γ*problem.T^α)
    uT = Χ*diagm(modified_spectrum)*(Χ/problem.u0)
    return MittagLefflerEigenSolution(uT)
end