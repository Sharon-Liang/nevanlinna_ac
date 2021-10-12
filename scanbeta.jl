using Pkg; Pkg.activate("./")
using nevanlinna_ac
using DelimitedFiles
using Printf

println("2021-10-12: g=1.0, scan beta (10,20)")
beta = [i for i in range(10,20,step=0.1)]
len = length(beta)
chi8 = similar(beta)
chi28 = similar(beta)
η = 0.001

for i = 1:len
    β = beta[i]
    path1 = @sprintf "./data/scan_beta/g8_b_%.1f.txt" β
    path2 = @sprintf "./data/scan_beta/g28_b_%.1f.txt" β

    g8 = make_input(path1)
    chi8[i] = spectrum_density(0,η,g8)

    g28 = make_input(path2)
    chi28[i] = spectrum_density(0,η,g28)
end

open("./data/scan_beta/chi.txt","w") do file
    for i = 1:len
        writedlm(file,[beta[i] chi8[i] chi28[i]])
    end
end