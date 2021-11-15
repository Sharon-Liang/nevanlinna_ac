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


"""load Zi-Long Li data"""
g = 2.0
(d,r) = divrem(g, 1)
invT = [10, 20, 30, 40]
N = 40

plc = Vector{String}(undef, length(invT))
oplc = Vector{String}(undef, length(invT))
for i = 1:length(invT)
    plc[i] = @sprintf "./data/TFIsing-Li/imaginary_time_data/g_%ip%i_beta_%i_tau.txt" d 10*r invT[i]
    oplc[i] = @sprintf "./data/TFIsing-Li/imaginary_time_data/g_%ip%i_beta_%i_iwn.txt" d 10*r invT[i]
end

for i = 1:length(invT)
    open(oplc[i], "w") do file
        for n = 1:N
            d = readdlm(plc[i])
            res = giwn(n,invT[i],d[:,1],d[:,2])
            writedlm(file, [Masubara_freq(n,invT[i]) real(res) imag(res)])
        end
    end
end
        




"""setups"""
omega = [i for i in range(-4π,4π,step = π/400)]

"""Nevanlinna"""
pcollect = Vector{String}(undef, length(invT))
for i = 1:length(invT)
    pcollect[i] = @sprintf "./data/gt/gt28_g_%.1f_b_%i.txt" g invT[i]
end

pcollect2 = Vector{String}(undef, length(invT))
for i = 1:length(invT)
    pcollect2[i] = @sprintf "./data/gt/gt8_g_%.1f_b_%i.txt" g invT[i]
end

dc = Vector{Vector}(undef,length(invT))
dc2 = Vector{Vector}(undef,length(invT))
dl = Vector{Matrix}(undef,length(invT))
for i = 1:length(invT)
    dc[i] = chi_div_w(omega, pcollect[i], 5)
    dc2[i] = chi_div_w(omega, pcollect2[i], 5)
    dl[i] = chi_div_w(plca[i])
    d = readdlm(plcs[i]); bs0 = invT[i]*d[800,2]
    dl[i][800,2] = bs0
end

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

for i = 1:4
    s = sum(dl[i])* dω 
    str = @sprintf "β=%i, sum = %.5f" invT[i] s
    println(str)
end



β = 20; num = div(β, 10) |> Int

plot()
v = 2 * abs(g-1)
plot!([-v,v], seriestype="vline",line=(:dash, :black), label=false)
plot!(dl[num][:,1], dl[num][:,2], line=(:orange, 2),label="Li")
plot!(omega, dc[num],line=(:dash, 2),label="cmpo:n=5")
plot!(xlim=(-2π, 2π))
plot!(title = @sprintf "g = %.1f, β=20, η=0.001" g)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f_b20.pdf" g
savefig(spath)


nrange = [3,5,7,9]
li = [:dash, :solid]
c = palette(:tab20c)
dn = Vector{Vector}(undef,length(nrange))
for i = 1:length(nrange)
    dn[i] = chi_div_w(omega, pcollect[num], nrange[i])
end
plot!(dl[num][:,1], dl[num][:,2], line=(:black, 2),label="Li")
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




"""
de = [spectral_density(1.0, i, β) for i in omega]
plot!(omega, de ./omega, line=(:black,2),label="theory")
plot!(dl[num][:,1], dl[num][:,2], line=(:dash,:orange, 2),label="Li")
plot!(omega, dc[num],line=(:dash, 2),label="cmpo:n=5")
plot!(xlim=(-1,1))
#plot!(xlim=(0,2π), ylim=(0, maximum(de)*1.1))
plot!(title = @sprintf "g = %.1f, β = %i, η=0.001" g β)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f_b20.pdf" g
savefig(spath)

c = palette(:okabe_ito)
for j=1:4
    plot!(omega, dc[j], line=(1.7), label=@sprintf "cmpo:β=%i" invT[j])
end
for j=1:4
    plot!(dl[j][:,1], dl[j][:,2], line=(:dash, c[4+j],1), 
    marker=(:circle, 3, :white, stroke(1,c[4+j],1.0)), 
    label=@sprintf "Li:β=%i" invT[j])
end
plot!(xlim=(-0.5, 0.5))
plot!(title = @sprintf "g = %.1f, cmpo n=5, η=0.001" g)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f.pdf" g
savefig(spath)
"""

"""sum rule check"""





