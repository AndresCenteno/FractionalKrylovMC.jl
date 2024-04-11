"""
Checking in 1D whether the reparametrization trick for the exponential works
"""

using Statistics, QuadGK, Distributions, ForwardDiff

p(t,λ) = λ*exp(-t*λ)
A = sqrt(3)

Xdet(λ) = quadgk(t->exp(-t*A)*p(t,λ),0,Inf)[1]
Xdet(1.)
ForwardDiff.derivative(λ->Xdet(λ),1.)

Xsto(λ) = exp(-rand(Exponential(1/λ))*A)
mean([Xsto(1.) for _ in 1:1e4])

Xsto(λ,u) = exp(-(-log(1-u)/λ)*A)
expRNG(λ,u) = -log(1-u)/λ
mean([Xsto(1.,rand()) for _ in 1:1e4])

dXstodλ(λ,u) = exp(-(-log(1-u)/λ)*A)*(-A*log(1-u)/(λ^2))
mean([dXstodλ(1.,rand()) for _ in 1:1e4]) # esto esta funcionando

mean([ForwardDiff.derivative(λ->Xsto(λ),1.) for _ in 1:1e4])
# this works and this is our goal
mean([-A*Xsto(1.,u)*ForwardDiff.derivative(λ->expRNG(λ,u),1.) for u in rand(Int(1e4))])
#######################################################

"""
Checking in 1D whether this works for the AlphaStableSubordinator
"""

A = 3
Xdet(γ) = exp(-A^γ)
dXdetdγ2(γ) = exp(-A^γ)*log(A)*(-A^γ)
dXdetdγ(γ) = ForwardDiff.derivative(γ->Xdet(γ),γ)
Xdet(0.8)
dXdetdγ(0.8)
dXdetdγ2(0.8)

Xsto(γ,u,v) = exp(-A*stblrndsub(γ,u,v))
mean([Xsto(0.8,rand(),rand()) for _ = 1:1e6])
dXstodγ(γ,u,v) = -A*Xsto(γ,u,v)*ForwardDiff.derivative(γ->Xsto(γ,u,v),γ)
nsims = Int(1e6); U = rand(nsims,2)
mean([dXstodγ(0.8,U[i,1],U[i,2]) for i = 1:nsims])
Xsto(0.9,rand(),rand())
dXstodγ(0.8,rand(),rand())
ForwardDiff.derivative(γ->Xsto(γ,rand(),rand()),γ)
dXstodϵ(γ,u,v,ϵ) = (Xsto(γ+ϵ,u,v)-Xsto(γ-ϵ,u,v))/(2*ϵ); dXstodϵ(γ,u,v) = dXstodϵ(γ,u,v,1e-8)
mean([-A*Xsto(0.8,U[i,1],U[i,2])*dXstodϵ(0.8,U[i,1],U[i,2]) for i = 1:nsims])

# it should be working fuck, maybe quadrature?
using HCubature
exp(-A^0.8)
fracexpA(A,γ) = HCubature.hcubature(u->exp(-stblrndsub(γ,u[1],u[2])*A),[0;0],[1;1])[1]
ForwardDiff.derivative(γ->fracexpA(A,γ),0.8)

dfracexpAdγ(A,γ) = HCubature.hcubature(u->exp(-stblrndsub(γ,u[1],u[2])*A)*(-A)*ForwardDiff.derivative(γ->stblrndsub(γ,u[1],u[2]),γ),[0;0],[1;1])[1]
dfracexpAdγ(A,0.8)

γ = 0.8
dfracexpAdγMC = zeros(nsims)
for i=1:nsims
    u = rand(2)
    dfracexpAdγMC[i] = exp(-stblrndsub(γ,u[1],u[2])*A)*(-A)*ForwardDiff.derivative(γ->stblrndsub(γ,u[1],u[2]),γ)
end
mean(dfracexpAdγMC)
var(dfracexpAdγMC)