using nevanlinna_ac
using PyPlot
using DelimitedFiles
using Printf

fname = "delta"
path = @sprintf "./data/i%s.txt" fname
opath = @sprintf "./data/o%s_julia.txt" fname
g = make_input(path);

η = 0.001
omega = [i for i in range(-10,10,length=6000)];
A = [spectrum_density(w,η,g) for w in omega] 

open(opath, "w") do file
    for i = 1:6000
        writedlm(file, [omega[i] A[i]])
    end
end

"""
phis = core(g)
open("./data/phi_julia.txt","w") do file
    for i = 1:length(phis)
        writedlm(file, [i real(phis[i]) imag(phis[i])])
    end
end
"""