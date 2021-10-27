using Plots; pyplot(xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=12)
default(palette = palette(:okabe_ito))
using DelimitedFiles
using Printf
#using StatsFuns, SpecialFunctions
using nevanlinna_ac

"""
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
"""

function chi_div_w(w::Vector, path::String, n::Int64)
	η = 0.001
	g = make_input(path) 
	gs = MasubaraGF(n, g.GF[length(g.GF)-n+1:end])
	return [spectrum_density(i,η,gs) for i in w]
end

#setups
#o = [-i for i in range(-2π,2π,step = 4π/400)]
o = [-i for i in range(-1,1,step = 0.01)]
#brange = [10.0, 20.0, 30.0, 40.0]
#de = [spectral_density(1.0, i, 20) for i in o]

g = 1.0
β = 20

#load Zi-Long Li data
#pathl = @sprintf "./data/TFIsing-Li/Aw_%.1f_b20.txt" g
#dl = readdlm(pathl)





#n = 5 #number of Masubara frequencies used
nrange = [6, 7, 8]
chi = Array{Vector}(undef,3,2)
for  i = 1:3
    path1 = @sprintf "./data/gt/gt8_g_%.1f_b_%i.txt" g β
    path2 = @sprintf "./data/gt/gt28_g_%.1f_b_%i.txt" g β
    chi[i,1] = chi_div_w(o, path1, nrange[i])
    chi[i,2] = chi_div_w(o, path2, nrange[i])
end

#sc = 0.01 .* sum.(chi)
#for i = 1:4
#    s = @sprintf "%.5f" sc[i]
#    println(s)
#end


#v = 2 * abs(g-1)
plot()
#plot!(o, de ./ o, line=(:black,2),label="theory")
#plot!([-v,v], seriestype="vline",line=(:dash, :black), label=false)
#plot!(dl[:,1], dl[:,2] ./ dl[:,1], line=(:black,2), label="Li")
for j=1:3
    plot!(o, chi[j,1], line=(2), label=@sprintf "D=8,n=%i" nrange[j])
end
for j=1:3
    plot!(o, chi[j,2], line=(:dash, 1), 
    marker=(:circle, 2, stroke(0)), 
    label=@sprintf "D=2×8,n=%i" nrange[j])
end


plot!(xlim=(-1, 1))
plot!(title = @sprintf "g = %.1f, β=20, η = 0.001" g)
plot!(xlabel="ω", ylabel="χ''(ω)/ω")
spath = @sprintf "./notes/chi_g_%.1f_b_20_dchange.pdf" g
savefig(spath)
