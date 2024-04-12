using Revise
using MittagLeffler, ForwardDiff, FractionalKrylovMC, HCubature
"""
Full reparametrization trick with spectral kernel too
"""

A = 1/π
γ = 1/sqrt(2); α = 1/sqrt(3)
detres(α,γ,A) = mittleff(α,-A^γ)
graddetres(α,γ,A) = ForwardDiff.gradient(p->detres(p[1],p[2],A),[α;γ])

detres(α,γ,A)
detgrad = graddetres(α,γ,A)
n = ceil(γ/α)
μ = γ/(α*n)
Aμn = A^n

stores(α,γ,A) = hcubature(u->exp(-A^n*SpectralKernelRNG(α,1,u[1])^(α*n/γ)*stblrndsub(γ/(α*n),u[2],u[3])),[0;0;0],[1;1;1])[1]
stores(α,γ,A) # almost exact

TotalRNG(p,u) = SpectralKernelRNG(p[1],1,u[1])^(p[1]*n/p[2])*stblrndsub(p[2]/(p[1]*n),u[2],u[3])
dTotalRNGdp(p,u) = ForwardDiff.gradient(p->TotalRNG(p,u),p)
gradstores(α,γ,A) = hcubature(u->(-A^n)*dTotalRNGdp([α;γ],u)*exp(-A^n*SpectralKernelRNG(α,1,u[1])^(α*n/γ)*stblrndsub(γ/(α*n),u[2],u[3])),[0;0;0],[1;1;1])[1]
stograd = gradstores(α,γ,A)

@show detgrad, stograd