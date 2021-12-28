using nevanlinna_ac
using DelimitedFiles
using Printf
using Plots; gr(xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
default(palette = palette(:okabe_ito))

function tomatrix(mv::AbstractVector)
    M = length(mv)
    res = zeros(Ftype, M, 8)
    for i = 1:M, j = 1:4
        A = mv[i] |> transpose
        res[i,2j-1] = A[j] |> real
        res[i,2j] = A[j] |> imag
    end
    return res
end

function tocomplex(re, im)
end

# n in [5, 21, 40]
n = 21
#function ploss(preci::Number, n::Integer)
#    setprecision(preci)
    giwn_path = "./data/gaussian/giwn_gaussian_1.txt"
    sparam_path = @sprintf "./data/gaussian/sparam_gaussian_1_N_%i.txt" n
    abcd_path = @sprintf "./data/gaussian/abcd_gaussian_1_N_%i.txt" n
    theta_path = Vector{String}(undef, n-1)
    for i = 1:n-1
        theta_path[i] = @sprintf "./data/gaussian/theta_%i_gaussian_1_N_%i.txt" i n
    end
    s = readdlm(sparam_path)
    cs = (s[:,2] .+ one(Ctype)im .* s[:,3]) 
    abcd = readdlm(abcd_path)
    theta = [readdlm(theta_path[i]) for i=1:n-1]
    
    x1, y1 = readGF(giwn_path, num=n)
    x1n, y1n = toGeneralizedSchurdata(x1, y1, :f)
    println(eltype(x1n), " , ", eltype(y1n))
    ϕ, fac, prod = schur_parameter(x1n, y1n)

    fac = tomatrix(fac)
    prod = [tomatrix(prod[i,:]) for i=1:n-1]


    ds = ϕ .- cs

    dθ = similar(theta)
    for i = 1:n-1
        dθ[i] = similar(theta[i])
        dθ[i][:,1] = theta[i][:,1]
        dθ[i][:,2:end] = abs.(prod[i][i:end,:].- theta[i][:,2:end])
    end

    max_da = [maximum(abs.(dθ[i][:,2:3])) for i=1:n-1]
    max_db = [maximum(abs.(dθ[i][:,4:5])) for i=1:n-1]
    max_dc = [maximum(abs.(dθ[i][:,6:7])) for i=1:n-1]
    max_dd = [maximum(abs.(dθ[i][:,8:9])) for i=1:n-1]
#    return sp[:,1], abs.(dimg), abs.(fac_out .- abcd)
#end

x1, y1, dfac1 = ploss(128, 21)
x2, y2fac = ploss(2000, 21)

function plot_fac_diff(x1::AbstractVector, fac::AbstractMatrix)
    x1 = x1[1:end-1]
    plot(x1, fac[:,2], marker=(:circle, 3), line=(:dash,1), label="a.real")
    plot!(x1, fac[:,3], marker=(:circle, 3), line=(:dash,1), label="a.imag")
    plot!(x1, fac[:,4], marker=(:circle, 3), line=(:dash,1), label="b.real")
    plot!(x1, fac[:,4], marker=(:circle, 3), line=(:dash,1), label="a.real")
end
plot(x1[1:end-1], dfac1[:,2:3], marker=(:circle, 3), line=(:dash,1))


plot(x1, y1, lw = 2, label="128")
plot!(x2, y2, line=(:dash,2), label="2000")
plot!(yaxis=:log)

d1 = readdlm("./data/gaussian/sparam_gaussian_1_N_21.txt")
d2 = readdlm("./data/gaussian/sparam_gaussian_1_N_21_np.txt")

y = d1[:,2] .- d2[:,2] + 1im*(d1[:,3] - d2[:,3])
y = abs.(y)

plot(d1[3:end,1], y[3:end], yaxis=:log)

