module nevanlinna_ac
__precompile__()

using LinearAlgebra
using Random; Random.seed!()

#input file format: ωn  Re[G(iωn)]  Im[G(iωn)]
include("utilities.jl")

#function nevanlinna(path::String, num::Int64; type::Symbol = :f)
    data = readdlm(path)
    fnum = min(num, length(data[:,1])) # number of frenquencies used
    Y = [1.0im * data[i,1] for i=1:fnum] # iωn
    G = [data[i,2]+1.0im*data[i,3] for i=1:fnum] # G(iωn)

    check = isNevanlinnasolvable(Y, G)
    println("isNevanlinnasolvable: ", check[1], "; err = ", check[2])
    NG = -G

    
 #end



end
