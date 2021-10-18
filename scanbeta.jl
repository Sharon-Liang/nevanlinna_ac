using Pkg; Pkg.activate("./")
using nevanlinna_ac
using DelimitedFiles
using Printf

println("2021-10-14: scan g and beta")
η = 0.001
beta = [10.0, 20.0, 30.0, 40.0]
gamma = [0.8, 0.9, 1.0,  1.1, 1.2]
omega = [i for i in range(-5,5,length=1000)]

for j = 1:length(gamma), b = 1:length(beta)
    g = gamma[j]; β = beta[b]
    path1 = @sprintf "./data/gt/gt8_g_%.1f_b_%i.txt" g β
    path2 = @sprintf "./data/gt/gt28_g_%.1f_b_%i.txt" g β
    opath1 = @sprintf "./data/gt/a_g_%.1f_b_%i.txt" g β
    opath2 = @sprintf "./data/gt/a5_g_%.1f_b_%i.txt" g β

    g8 = make_input(path1); g8s = make_input(path1, 5)
    g28 = make_input(path2); g28s = make_input(path2, 5)

    chi8 = [spectrum_density(i,η,g8) for i in omega]
    chi28 = [spectrum_density(i,η,g28) for i in omega]

    chi8s = [spectrum_density(i,η,g8s) for i in omega]
    chi28s = [spectrum_density(i,η,g28s) for i in omega]

    open(opath1,"w") do file
        for i = 1:length(omega)
            writedlm(file,[omega[i] chi8[i] chi28[i]])
        end
    end

    open(opath2,"w") do file
        for i = 1:length(omega)
            writedlm(file,[omega[i] chi8s[i] chi28s[i]])
        end
    end
end