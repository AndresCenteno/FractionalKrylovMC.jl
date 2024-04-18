# checking where the NaNs come from

using FractionalKrylovMC, HCubature

flag = 1
k = 30

## KRYLOV
while true
    @show flag
    flag += 1
    global problem = create_random_problem(rand()*0.9+0.1,rand()*0.9+0.1,k,MittagLefflerRand())
    try
        sol = solve(problem,10,QuadKrySolver(),rtol=1e-2)
    catch
        break
    end
end
sol = solve(problem,QuadKrySolver())
# STANDARD
while true
    @show flag
    flag += 1
    global problem = create_random_problem(rand()*0.9+0.1,rand()*0.9+0.1,k,MittagLefflerRand())
    sol = solve(problem,QuadSolver(),rightΔ=1e-5)
    if any(isnan.(sol.uT)); break; end
end


sol = solve(problem,QuadSolver(),rightΔ=1e-5)
sol.uT
rightΔ = 1e-5
Anμ = problem.Anμ; n = problem.nμ; γ = problem.γ; α = problem.α; t = problem.t; u0 = problem.u0
TotalRNG(p,u) = SpectralKernelRNG(p[1],t^p[1],u[1])^(p[1]*n/p[2])*stblrndsub(p[2]/(p[1]*n),u[2],u[3])
dTotalRNGdp(p,u) = ForwardDiff.gradient(p->TotalRNG(p,u),p)
buffer = hcubature_buffer(u->exp(-Anμ*TotalRNG([α,γ],u))*u0,[0;0;0],(1-rightΔ).*[1;1;1])
hcubature(u->exp(-Anμ*TotalRNG([α,γ],u))*u0,[0;0;0],(1-rightΔ).*[1;1;1],rtol=1e-2,atol=1e-5,buffer=buffer)[1]

using DataStructures
collect_buffer = extract_all!(buffer)
typeof(collect_buffer[1])
fieldnames(HCubature.Box)
collect_buffer[1].kdiv

using Random
Random.seed!(1234)
A = randn(30,30); A = A*A'
T = prevfloat(Inf)
exp(-A*T)