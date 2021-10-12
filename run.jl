using nevanlinna_ac
using PyPlot
using DelimitedFiles
using Printf

fac = [0, 1]
A = x-> x*gaussian(x,fac)
x = [i for i in range(-5,5,length=1000)]
y = [A(i) for i in x ]

plot(x,y)
PyPlot.display_figs()

open("./data/ag3_bose_b20.txt","w") do file
    for i = 1:1000
        writedlm(file,[x[i] A(x[i])])
    end
end

num = 30
β = 20


B = x -> gaussian(x,fac)
ωn = [Masubara_freq(i,β,type=:b) for i=1:num]
GF = [Masubara_GF(i,B,β,type=:b) for i=1:num]

open("./data/ig3_bose_b20.txt","w") do file
    for i = 1:num
        writedlm(file,[ωn[i] real(GF[i]) imag(GF[i])])
    end
end