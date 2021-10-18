using PyPlot
using DelimitedFiles
using Printf

g = 0.9
beta = [10.0, 20.0, 30.0, 40.0]

figure()
for β in beta
    opath2 = @sprintf "./data/gt/a5_g_%.1f_b_%i.txt" g β
    d = readdlm(opath2)
    llab = @sprintf "β = %i" β
    plot(d[:,1], d[:,3], lw = 2, label=llab)
end
xlim(-0.1,0.1)
xlabel("ω")
ylabel("χ(ω)/ω")
title(@sprintf "g = %.1f, N = 5" g)
fname = @sprintf "./g_%.1f_n_5.pdf" g
savefig(fname)
