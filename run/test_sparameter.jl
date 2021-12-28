using nevanlinna_ac
using DelimitedFiles
using Printf
using Plots; gr(xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
default(palette = palette(:okabe_ito))

#setprecision(256)
# n in [5, 21, 40]
function ploss(preci::Number, n::Integer)
    setprecision(preci)
    path = "./data/gaussian/giwn_gaussian_1.txt"
    spath = @sprintf "./data/gaussian/sparam_gaussian_1_N_%i.txt" n
    sp = readdlm(spath)
    x1, y1 = readGF(path, num=n)
    x1n, y1n = toGeneralizedSchurdata(x1, y1, :f)

    # to generalized schur data
    println(eltype(x1n), " , ", eltype(y1n))
    ϕ = schur_parameter(x1n, y1n)

    dimg = ϕ .- (sp[:,2] .+ one(Complex{BigFloat})im .* sp[:,3]) 
    return sp[:,1], abs.(dimg)
end


x1, y1 = ploss(128, 21)
x2, y2 = ploss(2000, 21)
 
plot(x1, y1, lw = 2, label="128")
plot!(x2, y2, line=(:dash,2), label="2000")
plot!(yaxis=:log)

d1 = readdlm("./data/gaussian/sparam_gaussian_1_N_21.txt")
d2 = readdlm("./data/gaussian/sparam_gaussian_1_N_21_np.txt")

y = d1[:,2] .- d2[:,2] + 1im*(d1[:,3] - d2[:,3])
y = abs.(y)

plot(d1[3:end,1], y[3:end], yaxis=:log)

