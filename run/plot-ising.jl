using Plots; pyplot(xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=12)
default(palette = palette(:okabe_ito))
using DelimitedFiles
using Printf
using StatsFuns, SpecialFunctions
using nevanlinna_ac

"""
    Check G(iω_n)
"""
function Masubara_freq(n::Int64, β::Real; type::Symbol=:b)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end

function giwn(n::Int, β::Real, τ::Vector, g::Vector)
    isapprox(β, τ[end], atol=0.1) ? iωn = 1.0im*Masubara_freq(n,β,type=:b) : error("wrong β")
    len = length(τ)
    res = 0.0 + 0.0im
    for i = 1: len
        res += -g[i] * exp(iωn*τ[i])
    end
    return β/len*res
end

function spectral_density(J::Real, ω::Real, β::Real; η::Float64=0.05,Γ::Real=1.)
    # 2Imχ(ω)
    if J/Γ == 0.
        return 2π*sinh(Γ*β)/cosh(Γ*β)*(delta(ω-2Γ,η) - delta(ω+2Γ,η))
    elseif J/Γ == 1.
        # Ref: PRX.4.031008(2014) A6
        g0 = 0.858714569; zc = J^(-1/4)
        T = 1/β
        up = zc * g0 * β^(3/4)
        down = 2^(1/4) * √π * gamma(1/8) * gamma(5/8)
        fac = up/down
        res = sinh(ω/(2*T)) * abs(gamma(1/8 - 1im*ω/(2π*T)))^2
        return 2*fac*res
    else
        @error "No exact results. J/Γ should be 1. or 0."
    end
end

function chi_div_w(path::String)
    dl = readdlm(path)
    dl[:,2] = dl[:,2] ./ dl[:,1]
    return dl
end

function chi_div_w(w::Vector, path::String, n::Int64)
	η = 0.001
	g = make_input(path) 
	gs = MasubaraGF(n, g.GF[length(g.GF)-n+1:end])
	return [spectrum_density(i,η,gs) for i in w]
end


"""Nevanlinna"""
omega = [i for i in range(-4π,4π,step = π/400)]
#load Zi-Long Li data
g = 1.0
(d,r) = divrem(g, 1)
invT = [10, 20, 30, 40]
N = 40

plca = Vector{String}(undef, length(invT))
plcs = Vector{String}(undef, length(invT))
for i = 1:length(invT)
    plca[i] = @sprintf "./data/ising-Li/g_%ip%i_beta_%i_A.txt" d 10*r invT[i]
    plcs[i] = @sprintf "./data/ising-Li/g_%ip%i_beta_%i_S.txt" d 10*r invT[i]
end

#load ls data
D = [8, 16]
pc = Matrix{String}(undef, length(invT), 2*length(D))
for i = 1:length(invT), j=1:length(D)
    pc[i,2j-1] = @sprintf "./data/ising-imagtime/gdivwn_g_%.1f_D_%i_beta_%i.txt" g D[j] invT[i]
    pc[i,2j] = @sprintf "./data/ising-imagtime/gdivwn_g_%.1f_D_2m%i_beta_%i.txt" g D[j] invT[i]
end

dc = Matrix{Vector}(undef, length(invT), 2*length(D))
dl = Vector{Matrix}(undef,length(invT))

n = 5
for i = 1:length(invT)
    dl[i] = chi_div_w(plca[i])
    d = readdlm(plcs[i]); bs0 = invT[i]*d[800,2]
    dl[i][800,2] = bs0
    for j=1:2*length(D)
        dc[i,j] = chi_div_w(omega, pc[i,j], n)
    end
end

for i = 1:length(invT)
    path = @sprintf "./data/ising-Li/g_%ip%i_beta_%i_Adivw.txt" d 10*r invT[i]
    open(path,"w") do file
        for w = 1:size(dl[i])[1]
            writedlm(file, [dl[i][w,1] dl[i][w,2]])
        end
    end
end

β = 20
t = div(β,10)

opc = Matrix{String}(undef, length(invT), 2*length(D))
for i = 1:length(invT), j=1:length(D)
    opc[i,2j-1] = @sprintf "./data/ising-spectrum/Aw_ac_g_%.1f_D_%i_beta_%i_n_%i.txt" g D[j] invT[i] n
    opc[i,2j] = @sprintf "./data/ising-spectrum/Aw_ac_g_%.1f_D_2m%i_beta_%i_n_%i.txt" g D[j] invT[i] n
end

for i = 1:length(invT), j=1:2*length(D)
    open(opc[i,j],"w") do file
        for w=1:length(omega)
            out = dc[i,j] .* omega
            writedlm(file, [omega[w] out[w]])
        end
    end
end




#g!=1.0
plot()
v = 2 * abs(g-1)
plot!([-v,v], seriestype="vline",line=(:dash, :black), label=false)
plot!(dl[t][:,1], dl[t][:,2], line=(:solid, 2),label="Li")
plot!(omega, dc[t,1],line=(:dash, 2),label="χ=8")
plot!(omega, dc[t,2],line=(:dash, 2),label="χ=8×2")
plot!(omega, dc[t,3],line=(:dash, 2),label="χ=16")
plot!(omega, dc[t,4],line=(:dash, 2),label="χ=16×2")
plot!(xlim=(-2π, 2π))
#plot!(xlim=(0,2π), ylim=(0, maximum(de)*1.1))
plot!(title = @sprintf "g=%.1f, β=%i, η=0.001, n=%i" g β n)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f_b%i_n%i.pdf" g invT[t] n
savefig(spath)


#g=1.0
de = [spectral_density(1.0, i, β) for i in omega]
for t = 1:length(invT)
    de = [spectral_density(1.0, i, invT[t]) for i in omega]
    path = @sprintf "./data/ising-spectrum/Aw_exact_g_%.1f_beta_%i.txt" g invT[t]
    open(path,"w") do file
        for w = 1:length(omega)
            writedlm(file, [omega[w] de[w]])
        end
    end
end

plot(omega, de ./omega, line=(:black,2),label="theory")
plot!(dl[t][:,1], dl[t][:,2], line=(:solid, 2),label="Li")
plot!(omega, dc[t,1],line=(:dash, 2),label="χ=8")
plot!(omega, dc[t,2],line=(:dash, 2),label="χ=8×2")
plot!(omega, dc[t,3],line=(:dash, 2),label="χ=16")
plot!(omega, dc[t,4],line=(:dash, 2),label="χ=16×2")
plot!(xlim=(-0.25, 0.25))
#plot!(xlim=(0,2π), ylim=(0, maximum(de)*1.1))
plot!(title = @sprintf "g=%.1f, β=%i, η=0.001, n=%i" g β n)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f_b%i_n%i.pdf" g invT[t] n
savefig(spath)


#nchange
nrange = [3,5,7,9]
li = [:dash, :solid]
c = palette(:tab20c)
dn = Vector{Vector}(undef,length(nrange))
for i = 1:length(nrange)
    dn[i] = chi_div_w(omega, pc[t,], nrange[i])
end
plot!(dl[t][:,1], dl[t][:,2], line=(:black, 2),label="Li")
for j=1:length(nrange)
    plot!(omega, dn[j], line=(li[rem(j,2)+1], c[nrange[j]],1.5), 
    label=@sprintf "n=%i" nrange[j])
end 
plot!(xlim=(-1,1))
plot!(title = @sprintf "g = %.1f, β=20, η=0.001" g)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f_b20_nchange.pdf" g
savefig(spath)




nrange = [4,5,6,7]
dn = Vector{Vector}(undef,length(nrange))
dn2 = Vector{Vector}(undef,length(nrange))
for i = 1:length(nrange)
    dn[i] = chi_div_w(omega, pcollect[num], nrange[i])
    dn2[i] = chi_div_w(omega, pcollect2[num], nrange[i])
end
c = palette(:tab20c)
plot!(dl[num][:,1], dl[num][:,2], line=(:black, 2),label="Li")
for j=1:length(nrange)
    plot!(omega, dn[j], line=(:solid, c[4j-3], 2), 
    label=@sprintf "D=16, n=%i" nrange[j])
    plot!(omega, dn2[j], line=(:dash, c[4j-2],1.5), 
    marker=(:circle, 4, :white, stroke(1.5,c[4j-2],1.0)), 
    label=@sprintf "D=8 , n=%i" nrange[j])
end
plot!(xlim=(-0.5,0.5))
plot!(title = @sprintf "g = %.1f, β=20, η=0.001" g)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f_b20_dchange.pdf" g
savefig(spath)






"""sum rule check"""
dω = π/400
for i = 1:4
    s1 = sum(dc[i])* dω 
    s2 = sum(dc2[i])* dω
    str1 = @sprintf "β=%i, D=16, sum = %.5f" invT[i] s1
    str2 = @sprintf "β=%i, D=8, sum = %.5f" invT[i] s2
    println(str1)
    println(str2)
end



"""Check integral data"""
inputs = readdlm("./data/imagtime/gtau_g_1.0_D_8_beta_20.txt")
outputs = [giwn(i,20,inputs[:,1].*20,inputs[:,2]) for i = 1:40]
anly = readdlm("./data/imagtime/giwn_g_1.0_D_8_beta_20.txt")

c = palette(:okabe_ito)
scatter(anly[:,1], anly[:,2] .- real.(outputs), line=(:dash, c[5], 1), 
    marker=(:circle, 6, c[5], stroke(1.5,c[5],1.0)), 
    label="analytic")
scatter!(anly[:,1], real.(outputs), 
    marker=(:circle, 6, :white, stroke(1.5,c[6],1.0)), 
    label="integral")
plot!(xlabel="ωn", ylabel="Re ΔG(iωn)", title="g=1.0, β=20, χ=8",
    legend=:bottomright)
savefig("./notes/check_integral_real_diff.pdf")



"""load Zi-Long Li data"""
g = 1.0
(d,r) = divrem(g, 1)
invT = [10, 20, 30, 40]
N = 40

plc = Vector{String}(undef, length(invT))
oplc = Vector{String}(undef, length(invT))
for i = 1:length(invT)
    plc[i] = @sprintf "./data/TFIsing-Li/imaginary_time_data/g_%ip%i_beta_%i_tau.txt" d 10*r invT[i]
    oplc[i] = @sprintf "./data/TFIsing-Li/imaginary_time_data/g_%ip%i_beta_%i_iwn.txt" d 10*r invT[i]
end

"""load ls data"""
D = [8, 16]
pc = Matrix{String}(undef, length(invT), 2*length(D))
opc = Matrix{String}(undef, length(invT), 2*length(D))
for i = 1:length(invT), j=1:length(D)
    pc[i,2j-1] = @sprintf "./data/imagtime/gtau_g_%.1f_D_%i_beta_%i.txt" g D[j] invT[i]
    pc[i,2j] = @sprintf "./data/imagtime/gtau_g_%.1f_D_2m%i_beta_%i.txt" g D[j] invT[i]
    opc[i,2j-1] = @sprintf "./data/imagtime/giwn_g_%.1f_D_%i_beta_%i.txt" g D[j] invT[i]
    opc[i,2j] = @sprintf "./data/imagtime/giwn_g_%.1f_D_2m%i_beta_%i.txt" g D[j] invT[i]
end

"""compare G(iωn)"""
β = 40
t=div(β,10)

dl = readdlm(oplc[t])
dc = [readdlm(opc[t,i]) for i=1:2*length(D)]

c = palette(:okabe_ito)
#scatter(dl[:,1], dl[:,2] , line=(:dash, :black, 1), 
#    marker=(:circle, 6, :black, stroke(1.5,:black,1.0)), 
#    label="Li")
scatter(dc[1][:,1], dc[1][:,2] .- dl[:,2], 
    marker=(:circle, 6, stroke(0.5,1.0)), 
    label="χ=8")
scatter!(dc[2][:,1], dc[2][:,2] .- dl[:,2], 
    marker=(:circle, 6, stroke(0.5,1.0)), 
    label="χ=8×2")
scatter!(dc[3][:,1], dc[3][:,2] .- dl[:,2], 
    marker=(:circle, 6, stroke(0.5,1.0)), 
    label="χ=16")
scatter!(dc[4][:,1], dc[4][:,2] .- dl[:,2], 
    marker=(:circle, 6, stroke(0.5,1.0)), 
    label="χ=16×2")
plot!(xlabel="ωn", ylabel="|Re(ΔG(iωn))|", title="g=1.0, β=20, χ=8",
    legend=:bottomright)
savefig("./notes/giwn_real_diff.pdf")




"""compare G(τ)"""
β = 40
t=div(β,10)

dl = readdlm(plc[t])
dc = [readdlm(pc[t,i]) for i=1:2*length(D)]

for i = 1:2*length(D)
    dc[i][1:1600,2] = abs.(dc[i][1:1600,2] .- dl[:,2])
end

plot(dc[1][1:1600,1],dc[1][1:1600,2], line=(2),label="χ=8")
plot!(dc[2][1:1600,1],dc[2][1:1600,2], line=(2),label="χ=2×8")
plot!(dc[3][1:1600,1],dc[3][1:1600,2], line=(2),label="χ=16")
plot!(dc[4][1:1600,1],dc[4][1:1600,2], line=(2),label="χ=2×16")
ptitle = @sprintf "G(τ) difference of g=%.1f, at β=%i" g invT[t]
plot!(xlabel="τ/β", ylabel="|ΔG(τ)|", title=ptitle)
spath = @sprintf "./notes/gtau_diff_g_%.1f_beta_%i.pdf" g invT[t]
savefig(spath)

  











"""sum rule check"""





