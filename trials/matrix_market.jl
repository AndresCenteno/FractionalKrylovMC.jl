using MatrixMarket, LinearAlgebra, Arpack, FractionalKrylovMC

A = mmread("trials/50.mtx")
eigen(Matrix(A))
Adense = Matrix(A)
maximum(Adense)

B = mmread("trials/lapl27000_3D.mtx")

