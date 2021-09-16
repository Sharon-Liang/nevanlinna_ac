using DelimitedFiles
using LinearAlgebra
#delta funct"ion data
#d = readdlm("/Users/liangshuang/Desktop/CAS-code/nevanlinna_ac/References/run/idelta.txt")
d = readdlm("/Users/liangshuang/Desktop/CAS-code/nevanlinna_ac/References/run/ig3.txt")
fnum = length(d[:,1])
Y = [1.0im * d[i,1] for i=1:fnum]
G = [d[i,2]+1.0im*d[i,3] for i=1:fnum]

p = pick_matrix(Y, G)
check = isNevanlinnasolvable(Y, G)
println("isNevanlinnasolvable: ", check[1], "; minimum eval= ", check[2])

evals = eigvals(p)