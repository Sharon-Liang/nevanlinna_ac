using nevanlinna_ac
using DelimitedFiles
using Printf
using Plots; pyplot(xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
default(palette = palette(:okabe_ito))

setprecision(128)

function mti(z::Number) 
    return mt(z, 1.0im)
end

function imti(z::Number)
    return imt(z, 1.0im)
end

# n2 in [5, 21, 40]
n2 = 40
path = "./data/gaussian/giwn_gaussian_1.txt"
spath = @sprintf "./data/gaussian/sparam_gaussian_1_N_%i.txt" n2
sp = readdlm(spath)
x1, y1 = readGF(path, num=n2)
x1n, y1n = toNevanlinnadata(x1, y1, :f)

# to generalized schur data
y1n = mti.(y1n)
ϕ = schur_parameter(x1n, y1n)

plot(sp[:,1], sp[:,2], 
    line=(:dash, 1), marker=(:circle, 5, stroke(0.3)),
    label="c++")
plot!(imag.(x1n), real.(ϕ), 
    line=(:dash, 1), marker=(:star, 5, stroke(0.3)),
    label="julia")
plot!(xlabel="ωn", ylabel="imag Schur parameter")
f1 = @sprintf "./imag_sp_N_%i.pdf" n2
savefig(f1)


plot(sp[:,1], sp[:,3], 
    line=(:dash, 1), marker=(:circle, 5, stroke(0.3)),
    label="c++")
plot!(imag.(x1n), imag.(ϕ), 
    line=(:dash, 1), marker=(:star,5, stroke(0.3)),
    label="julia")
plot!(xlabel="ωn", ylabel="real Schur parameter")
f2 = @sprintf "./real_sp_N_%i.pdf" n2
savefig(f2)