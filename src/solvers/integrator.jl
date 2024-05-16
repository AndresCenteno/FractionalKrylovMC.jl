struct IntegratorSolver end

function solve(problem::MittagLefflerProblem{T},::IntegratorSolver; Nt = 100) where T
    Δt = 1/Nt
    list_of_times = zeros(Float64,Nt,Nt,Nt)
    Anμ = problem.Anμ; n = problem.nμ; γ = problem.γ; α = problem.α; t = problem.t; u0 = problem.u0
    TotalRNG(p,u) = SpectralKernelRNG(p[1],t^p[1],u[1])^(p[1]*n/p[2])*stblrndsub(p[2]/(p[1]*n),u[2],u[3])
    for i=1:Nt
        for j=1:Nt
            for k=1:Nt
                list_of_times[i,j,k] = TotalRNG([α;γ],[(i-1)*Δt;(j-1)*Δt;(k-1)*Δt])
            end
        end
    end
    list_of_times = sort(list_of_times[:])
    uT = zeros(length(u0))
    heat_T = u0
    for l=1:(Nt^3-1)
        if list_of_times[l] < 6
        heat_T = (I-Anμ*(list_of_times[l+1]-list_of_times[l]))\heat_T
        uT += heat_T
        end
    end
    return uT, list_of_times
end