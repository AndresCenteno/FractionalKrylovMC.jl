abstract type MatVecSolver end
abstract type MatVecSolution end
#TODO: perhaps merge MittagLeffler and just FracExp into one type

################################
### MITTAG-LEFFLER + FRAC ######
################################

"""
MittagLefflerSolution stores solution u(t) = E_α(-t^α*A^γ)*u0
and du(t)/dα, du(t)/dγ in a structure for convenience
"""
struct MittagLefflerSolution{T} <: MatVecSolution
    uT::Vector{T}
    duTdα::Vector{T}
    duTdγ::Vector{T}
    function MittagLefflerSolution(uT::Vector{T},duTdα::Vector{T},duTdγ::Vector{T}) where T<:Real
        new{T}(uT,duTdα,duTdγ)
    end
end

struct MittagLefflerProblem{T}
    N::Int # size of matrix
    A::Matrix{T}
    Anμ::Matrix{T}
    u0::Vector{T}
    t::T
    α::T # exponent in time
    #TODO
    # find if the spectral Kernel for β≠1 is a pdf, postponed indefinetly
    # β::T
    γ::T # exponent in space
    μ::T
    nμ::Int
    #TODO
    # include Krylov decomposition
    
    # and include spectral decomposition why not
    λs::Vector{T}
    xs::Matrix{T}
    solution::MittagLefflerSolution{T}

    function MittagLefflerProblem(A::Matrix{T},u0::Vector{T},t::T,α::T,γ::T) where T<:Real
        N = length(u0)
        if size(A,1) != size(A,2)
            error("Matrix provided must be squared")
        elseif α > 1 || α < 0
            error("α time exponent must be between 0 and 1")
        elseif γ > 1 || γ < 0
            error("γ space exponent must be between 0 and 1")
        elseif t <= 0
            error("time must be positive")
        elseif N != size(A,1)
            error("Matrix must have same size as vector")
        end
        nμ = Int(ceil(γ/α))
        μ = γ/(α*nμ)
        new{T}(N,A,A^nμ,u0,t,α,γ,μ,nμ,nothing,nothing,nothing)
    end

    function MittagLefflerProblem(λs::Vector{T},xs::Matrix{T},u0::Vector{T},t::T,α::T,γ::T) where T<:Real
        N = length(u0)
        A = xs*diagm(λs)*inv(xs) # ill posed xs*diagm(λs)/xs
        if size(A,1) != size(A,2)
            error("Matrix provided must be squared")
        elseif α > 1 || α < 0
            error("α time exponent must be between 0 and 1")
        elseif γ > 1 || γ < 0
            error("γ space exponent must be between 0 and 1")
        elseif t <= 0
            error("time must be positive")
        elseif N != size(A,1)
            error("Matrix must have same size as vector")
        end
        nμ = Int(ceil(γ/α))
        μ = γ/(α*nμ)
        uT_spectral = xs*diagm(mittleff.(α,-λs.^γ*t^α))*(xs\u0)
        ϵ = sqrt(eps())
        duTdα_FD = xs*diagm((mittleff.(α+ϵ,-λs.^γ*t^(α+ϵ))-mittleff.(α-ϵ,-λs.^γ*t^(α-ϵ)))/(2*ϵ))*(xs\u0)
        duTdγ_FD = xs*diagm((mittleff.(α,-λs.^(γ+ϵ)*t^(α))-mittleff.(α,-λs.^(γ-ϵ)*t^(α)))/(2*ϵ))*(xs\u0)
        new{T}(N,A,A^nμ,u0,t,α,γ,μ,nμ,λs,xs,MittagLefflerSolution(uT_spectral,duTdα_FD,duTdγ_FD))
    end
end

struct MittagLefflerMCSolution{T} <: MatVecSolution 
    uT::Union{Vector{T},Matrix{T}} # 
    times::Vector{T}
    nsims::Int
    cutoff::T
    # to update to store Σ and checks
    function MittagLefflerMCSolution(uT::Union{Vector{T},Matrix{T}},times::Vector{T},
                                    nsims::Int,cutoff::T) where T<:Real
        new{T}(uT,times,nsims,cutoff)
    end
end

struct MittagLefflerEigenSolution{T} <: MatVecSolution
    uT::Vector{T}
    # can we have an estimate of error?
end

####################
## EXP + FRAC ######
####################

struct FracExpProblem{T}
    N::Int
    A::Matrix{T}
    u0::Vector{T}
    t::T
    γ::T
    # for exact solution
    Aγ::Matrix{T}
    λs::Vector{T}
    xs::Matrix{T}
    uT_spectral::Vector{T}
    duTdγ_spectral::Vector{T}
    # constructors
    # without spectral decomposition
    function FracExpProblem(A::Matrix{T},u0::Vector{T},t::T,γ::T) where T<:Real
        N = length(u0)
        (γ < 0.1 || γ >= 1) && throw(DomainError(γ,"Fractional exponent must be between 0.1 and 1"))
        new{T}(N,A,u0,t,γ,nothing,nothing,nothing,nothing,nothing)
    end
    # with spectral decomposition
    function FracExpProblem(λs::Vector{T},xs::Matrix{T},u0::Vector{T},t::T,γ::T) where T<:Real
        N = length(u0)
        A = xs*diagm(λs)/xs
        (γ < 0.1 || γ >= 1) && throw(DomainError(γ,"Fractional exponent must be between 0.1 and 1"))
        Aγ = xs*diagm(λs.^γ)/xs
        uT_spectral = xs*diagm(exp.(-t*λs.^γ))*(xs\u0)
        duTdγ_spectral = xs*diagm(-t*log.(λs).*λs.^γ.*exp.(-t*λs.^γ))*(xs\u0)
        new{T}(N,A,u0,t,γ,Aγ,λs,xs,uT_spectral,duTdγ_spectral)
    end
end

struct FracExpMCSolution{T} <: MatVecSolution 
    uT::Union{Vector{T},Matrix{T}}
    duTdγ::Vector{T}
    times::Vector{T}
    nsims::Int
    cutoff::T
    # to update to store Σ and checks
    function FracExpMCSolution(uT::Union{Vector{T},Matrix{T}},times::Vector{T},
                                    nsims::Int,cutoff::T) where T<:Real
        new{T}(uT,nothing,times,nsims,cutoff)
    end
    function FracExpMCSolution(uT::Union{Vector{T},Matrix{T}},duTdγ::Vector{T},times::Vector{T},
                                    nsims::Int,cutoff::T) where T<:Real
        new{T}(uT,duTdγ,times,nsims,cutoff)
    end
end