using LinearAlgebra, ForwardDiff, AlphaStableDistributions, ExponentialAction

"""
Auxiliary function RNG(U,V,α) for α-stable subordinator
"""

function StableSubordinatorRNG(U, V, α, time) # do the promotion automatically plsss
    β=1; scale=cos(pi*α/2)^(1/α); loc=0
    (U < 0 || U > 1) && throw(DomainError(V, "random numbers must come from a uniform distribution"))
    (V < 0 || V > 1) && throw(DomainError(V, "random numbers must come from a uniform distribution"))
    (α < 0.1 || α > 2) && throw(DomainError(α, "α must be in the range 0.1 to 2"))
    abs(β) > 1 && throw(DomainError(β, "β must be in the range -1 to 1"))
    ϕ = (U - 0.5) * π
    if α == one(T) && β == zero(T)
        return loc + scale * tan(ϕ)
    end
    w = -log(V)
    α == 2 && (return loc + 2*scale*sqrt(w)*sin(ϕ))
    β == zero(T) && (return loc + scale * ((cos((1-α)*ϕ) / w)^(one(T)/α - one(T)) * sin(α * ϕ) / cos(ϕ)^(one(T)/α)))
    cosϕ = cos(ϕ)
    if abs(α - one(T)) > 1e-8
        ζ = β * tan(π * α / 2)
        aϕ = α * ϕ
        a1ϕ = (one(T) - α) * ϕ
        return loc + scale * (( (sin(aϕ) + ζ * cos(aϕ))/cosϕ * ((cos(a1ϕ) + ζ*sin(a1ϕ))) / ((w*cosϕ)^((1-α)/α)) ))
    end
    bϕ = π/2 + β*ϕ
    x = 2/π * (bϕ * tan(ϕ) - β * log(π/2*w*cosϕ/bϕ))
    α == one(T) || (x += β * tan(π*α/2))

    return (loc + scale * x)*time^α
end

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
A = zeros(1,1); A[1] = 3; u0 = ones(1); uT = exp.(-T*A.^γ).*u0; k = 1
# random time and initial condition
u0 = rand(k); T = rand()*2
# random exponent for matrix
γ = 0.1 + rand()*0.9; γ = 0.8

# storing matrix to the power of γ
Aγ = Q*diagm(Λ.^γ)*Q'

uT = exp(-T*Aγ)*u0
duTdγ_spectrum = Q*diagm(-T.*log.(Λ).*Λ.^γ.*exp.(-T*Λ.^γ))*Q'*u0
# ϵ = sqrt(eps()); duTdγ_finite_differences = (Q*(diagm(exp.(-T*Λ.^(γ+ϵ))-exp.(-T*Λ.^(γ-ϵ))))*Q'*u0)/(2*ϵ)

###
begin
    uT_MC_simple = zeros(k)
    uT_MC_complex = zeros(k)
    τ = AlphaStable(γ,1,cos(γ*pi/2)^(1/γ),0)
    nsims = Int(1e6); invnsims = 1/nsims
    for _ in 1:nsims
        time = rand(τ); time2 = StableSubordinatorRNG(rand(), rand(), γ, T)
        if !isnan(time) && time < 1e5 && time2 < 1e5            # really (A*T^(1/γ))^(γ) = A^γ*T
            uT_MC_simple += expv()
            uT_MC_simple += exp(-A*time*T^(1/γ))*u0*invnsims
            uT_MC_complex += exp(-A*time2*T^(1/γ))*u0*invnsims
        end
    end
end

uT_MC_simple
uT
relative_error(uT,uT_MC_complex)

begin
    uT_MC = zeros(k); uT_MC_package = zeros(k); duTdγ_MC = zeros(k)
    nsims = Int(1e6); cutoff = 1e7
    Uvec = rand(nsims); Vvec = rand(nsims) # RNG stream stored
    times = StableSubordinatorRNG.(Uvec, Vvec, γ, T); times_package = sort(T*(1/γ)*rand(AlphaStable(γ,1,cos(γ*pi/2)^(1/γ),0),nsims))
    #TODO: not evaluate twice lmao
    dtimesdγ = similar(times)
    for i in eachindex(dtimesdγ)
        dtimesdγ[i] = ForwardDiff.derivative(γ->StableSubordinatorRNG(Uvec[i], Vvec[i], γ, T), γ)
    end
    # filter times
    times = filter(x->x<=cutoff,filter(!isnan,times));
    times_package = filter(x->x<=cutoff,filter(!isnan,times_package))
    # sort times
    perm = sortperm(times); times = times[perm]; dtimesdγ = dtimesdγ[perm]
    jumps = [times[1];times[2:end].-times[1:end-1]]
    jumps_package = [times_package[1];times_package[2:end].-times_package[1:end-1]]
    v = copy(u0); v_package = copy(u0)
    for sim in eachindex(jumps)
        v = expv(jumps[sim],-A,v)
        uT_MC += v/nsims
        duTdγ_MC += dtimesdγ[sim].*v/nsims
    end
    for sim in eachindex(jumps_package)
        v_package = expv(jumps_package[sim],-A,v_package)
        uT_MC_package += v_package/nsims
    end
    duTdγ_MC = -A*duTdγ_MC
end

@show uT_MC_package
@show uT_MC
@show uT

@show duTdγ_spectrum
@show duTdγ_MC
