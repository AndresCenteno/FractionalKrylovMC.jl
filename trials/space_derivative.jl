using FractionalKrylovMC, ForwardDiff, LinearAlgebra, ExponentialAction

"""
Checking whether d(exp(-t*A^γ)*u0)/dγ = -t*log(A)*exp(-t*A^γ)*u0
can be stochastically computed throught the reparametrization trick as
exp(-t*A^γ)*u0=∫_0^1 exp(-t^{1/α}*p(α,u)*A)*u0 du so derivative should be
d(exp(-t*A^γ)*u0)/dγ = ∫_0^1 d(-t^{1/α}*p(α,u)*A)/dα * exp(-t^{1/α}*p(α,u)*A)*u0 du
"""
# size of matrix
k = 20 

# creating random spectral decomposition
Q = Matrix(qr(rand(k,k)).Q); Λ = rand(k)
A = Q*diagm(Λ)*Q'
# random time and initial condition
u0 = rand(k); Tf = rand()*2
# random exponent for matrix
γ = 0.1 + rand()*0.9

# storing matrix to the power of γ
Aγ = Q*diagm(Λ.^γ)*Q'

uT = exp(-Tf*Aγ)*u0
duTdγ_spectrum = Q*diagm(-Tf.*log.(Λ).*Λ.^γ.*exp.(-Tf*Λ.^γ))*Q'*u0
# ϵ = sqrt(eps()); duTdγ_finite_differences = (Q*(diagm(exp.(-T*Λ.^(γ+ϵ))-exp.(-T*Λ.^(γ-ϵ))))*Q'*u0)/(2*ϵ)

###

begin
    uT_MC = zeros(k); duTdγ_MC = zeros(k)
    nsims = Int(1e4); cutoff = 1e7
    Uvec = rand(nsims); Vvec = rand(nsims) # RNG stream stored
    times = stblrndsub(γ,Uvec,Vvec)*Tf^(1/γ)
    #TODO: not evaluate twice lmao
    dtimesdγ = similar(times)
    for i in eachindex(dtimesdγ)
        dtimesdγ[i] = ForwardDiff.derivative(γ->stblrndsub(γ,Uvec[i],Vvec[i])*Tf^(1/γ), γ)
    end
    # filter times
    times = filter(x->x<=cutoff,filter(!isnan,times));
    # sort times
    perm = sortperm(times); times = times[perm]; dtimesdγ = dtimesdγ[perm]
    jumps = [times[1];times[2:end].-times[1:end-1]]
    v = copy(u0)
    for sim in eachindex(jumps)
        v = expv(jumps[sim],-A,v)
        uT_MC += v/nsims
        duTdγ_MC += dtimesdγ[sim].*v/nsims
    end
    duTdγ_MC = -A*duTdγ_MC
end

@show uT_MC
@show uT

@show duTdγ_spectrum
@show duTdγ_MC
