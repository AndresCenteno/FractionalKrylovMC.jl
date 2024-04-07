abstract type MittagLefflerSolver end
abstract type MittagLefflerSolution end

struct MittagLefflerProblem{T}
    N::Int # size of matrix
    A::Matrix{T}
    Anμ::Matrix{T}
    u0::Vector{T}
    t::T
    α::T # exponent in time
    #TODO
    # find if the spectral Kernel for β≠1 is a pdf
    # β::T
    γ::T # exponent in space
    μ::T
    nμ::Int
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
        new{T}(N,A,A^nμ,u0,t,α,γ,μ,nμ)
    end
end

struct MittagLefflerMCSolution{T} <: MittagLefflerSolution 
    uT::Union{Vector{T},Matrix{T}} # 
    times::Vector{T}
    nsims::Int
    cutoff::T
    # to update to store Σ and checks
    function MittagLefflerMCSolution(uT::Vector{T},times::Vector{T},
                                    nsims::Int,cutoff::T) where T<:Real
        new{T}(uT,times,nsims,cutoff)
    end
end

struct MittagLefflerEigenSolution{T} <: MittagLefflerSolution
    uT::Vector{T}
    # can we have an estimate of error?
end