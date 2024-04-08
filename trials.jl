# auxiliary file for me to understand how to use Krylov-based methods

using ArnoldiMethod

ArnoldiWorkspace(Float64,100,20)

@edit ArnoldiWorkspace(Float16,1,1)