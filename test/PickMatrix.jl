using nevanlinna_ac, Test
using Printf

pcollect = Vector{String}(undef, 2)
pcollect[1] = "./data/gaussian/giwn_gaussian_1.txt"
pcollect[2] = "./data/gaussian/giwn_gaussian_3.txt"

for p in pcollect, n = 1:40 
    x, y = readGF(p, num = n)
    xn, yn = toGeneralizedSchurdata(x, y, :f)
    ϕ = schur_parameter(xn, yn)
    a = all( abs.(ϕ) .≤ 1 )
    b, _= isGeneralizedSchursovable(xn, yn, tolerance = 0.)
    a && b ? continue : 
        println("")
        println(p, ", n = ", n)
        println("   pick matrix positive semidefinite: ", b, "; |ϕ|≤ 1: ", a)
end
