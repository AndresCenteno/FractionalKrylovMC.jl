# stblrnd generator of stable random variables and
# stblsubrnd generator of stable subordinators

"""# taken from https://github.com/markveillette/stbl/blob/master/stblrnd.m
# modified to allow control of the RNG stream
# inputs from Mamikon Gulian's final thesis (not publicly available)
# calls stblrnd(a2, skewness,scaleg, center)
# a2 = alpha/2; he uses alpha = sqrt(3)
# skewness = 1;
# scaleg = (cos(pi*a2/2))ˆ(1/a2); 
# center = 0;
"""
function stblrnd(alpha::T,beta::T,gamma::T,delta::T,
    U1::Union{Vector{T},T},U2::Union{Vector{T},T}) where T<:Real
    # %STBLRND alpha-stable random number generator.
    # % R = STBLRND(ALPHA,BETA,GAMMA,DELTA) draws a sample from the Levy 
    # % alpha-stable distribution with characteristic exponent ALPHA, 
    # % skewness BETA, scale parameter GAMMA and location parameter DELTA.
    # % ALPHA,BETA,GAMMA and DELTA must be scalars which fall in the following 
    # % ranges :
    # %    0 < ALPHA <= 2
    # %    -1 <= BETA <= 1  
    # %    0 < GAMMA < inf 
    # %    -inf < DELTA < inf
    # %
    # %
    # % R = STBLRND(ALPHA,BETA,GAMMA,DELTA,M,N,...) or 
    # % R = STBLRND(ALPHA,BETA,GAMMA,DELTA,[M,N,...]) returns an M-by-N-by-... 
    # % array.   
    # % 
    # %
    # % References:
    # % [1] J.M. Chambers, C.L. Mallows and B.W. Stuck (1976) 
    # %     "A Method for Simulating Stable Random Variables"  
    # %     JASA, Vol. 71, No. 354. pages 340-344  
    # %
    # % [2] Aleksander Weron and Rafal Weron (1995)
    # %     "Computer Simulation of Levy alpha-Stable Variables and Processes" 
    # %     Lec. Notes in Physics, 457, pages 379-392
    # %
    
    # % Check parameters
    if alpha <= 0 || alpha > 2
        throw(DomainError(alpha,"Alpha must be a scalar which lies in the interval (0,2]"))
    end
    if abs(beta) > 1
        throw(DomainError(beta,"Beta must be a scalar which lies in the interval [-1,1]"))
    end
    if gamma < 0
        throw(DomainError(gamma,"Gamma must be a non-negative scalar"))
    end
        
    # %---Generate sample----
    
    # % See if parameters reduce to a special case, if so be quick, if not 
    # % perform general algorithm
    
    if alpha == 2                  # % Gaussian distribution 
        r = sqrt(2) * randn(size(U1,1));
    
    elseif alpha==1 && beta == 0   # % Cauchy distribution
        r = tan.( pi/2 * (2*U1 .- 1) ); 
    
    elseif alpha == .5 && abs(beta) == 1 # % Levy distribution (a.k.a. Pearson V)
        r = beta ./ randn(size(U1,1)).^2;
    
    elseif beta == 0                # % Symmetric alpha-stable
        V = pi/2 * (2*U1 .- 1); 
        W = -log.(U2);          
        r = sin.(alpha * V) ./ ( cos.(V).^(1/alpha) ) .* 
            ( cos.( V.*(1-alpha) ) ./ W ).^( (1-alpha)/alpha ); 
    
    elseif alpha != 1                # % General case, alpha not 1
        V = pi/2 * (2*U1 .- 1); 
        W = - log.(U2);       
        C = beta * tan.(pi*alpha/2);
        B = atan.(C);
        S = (1 + C*C).^(1/(2*alpha));
        r = S * sin.( alpha*V .+ B ) ./ ( cos.(V) ).^(1/alpha) .*
           ( cos.( (1-alpha) * V .- B ) ./ W ).^((1-alpha)/alpha);
    
    else                             # % General case, alpha = 1
        V = pi/2 * (2*U1 .- 1); 
        W = - log(U2);          
        piover2 = pi/2;
        sclshftV =  piover2 + beta * V ; 
        r = 1/piover2 * ( sclshftV .* tan.(V) .- beta * log.( (piover2 * W .* cos(V) ) ./ sclshftV ) );      
              
    end
        
    # % Scale and shift
    if alpha != 1
       return gamma * r .+ delta;
    else
       return gamma * r .+ (2/pi) * beta * gamma * log.(gamma) .+ delta;  
    end    
end

stblrnd(alpha::T,beta::T,gamma::T,delta::T,N::Int) where T<:Real = stblrnd(alpha,beta,gamma,delta,rand(T,N),rand(T,N))
stblrnd(alpha::T,beta::T,gamma::T,delta::T) where T<:Real = stblrnd(alpha,beta,gamma,delta,1)

# subordinators
stblrndsub(alpha::T,U1::Union{Vector{T},T},U2::Union{Vector{T},T}) where T<:Real = stblrnd(alpha,one(eltype(alpha)),cos(alpha*pi/2)^(1/alpha),zero(eltype(alpha)),U1,U2)

stblrndsub(alpha::T,N::Int) where T<:Real = stblrndsub(alpha,rand(T,N),rand(T,N))
stblrndsub(alpha::T) where T<:Real = stblrndsub(alpha,1)[1]