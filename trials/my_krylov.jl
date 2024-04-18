using MatrixMarket

A = mmread("trials/784-convect.mtx")
KryDim = 10
A = Matrix(A)
u0 = rand(size(A,1))
# exp(A)u0 = norm(u0)*exp(A)*(u0/norm(u0))
u00 = u0/norm(u0)
V, H = my_arnoldi(A,KryDim,u00)
det = exp(A)*u0
e1 = zeros(KryDim); e1[1] = 1
kry = norm(u0)*V*exp(H)*e1
norm(kry-det)
function my_arnoldi(A::Matrix{T},KryDim::Int,v::Vector{T}) where T<:Real
    n = size(A,1);
    V = zeros(n,KryDim)
    H = zeros(KryDim,KryDim)
    
    V[:,1] = v
    for k=2:KryDim+1
        w = A*V[:,k-1]
        for j=1:k-1
            b = V[:,j]
            dotProd = (w'*b)
            w = w .- dotProd*b
            H[j,k-1] = dotProd
        end
        if k != KryDim + 1
            H[k,k-1] = norm(w)
            V[:,k] = w./norm(w)
        end
    end
    return V, H
end