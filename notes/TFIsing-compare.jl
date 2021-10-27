### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 6d565956-314f-11ec-0005-a12d4cf1fb59
begin
	using Pkg; Pkg.activate("../")
	using nevanlinna_ac
	using DelimitedFiles
	using Plots
	gr(size=(600,500), 
		xtickfontsize=13, 
		ytickfontsize=13, 
		xguidefontsize=16, 
		yguidefontsize=16, 
		legendfontsize=12, 
		dpi=100,
		grid=(:y, :gray, :solid, 1, 0.4));
	using Printf
	using StatsFuns, SpecialFunctions
	"Packages"
end

# ╔═╡ 93711dc9-3b99-4ce3-8b39-570b4a833f58
η = 0.001

# ╔═╡ 25c4b6ca-c100-488e-8f73-e514f001aac8
dω = 4π/799

# ╔═╡ 1ed9dc6a-2971-4c61-aa3c-b3f91f7cb34b
md"""
### $g = 1.0, β = 20$
"""

# ╔═╡ eaa15cb3-8958-4967-b3b1-51f509b83b8a
md"""
### $g = 1.5, β = 20$
"""

# ╔═╡ 8f5723b4-5871-46b9-a706-495b1f34de77
md"""
### $g = 2.0, β = 20$
"""

# ╔═╡ 424f1e78-bbb1-4d55-8c14-a64d2f61ab0f
md"""
## check sum rule
"""

# ╔═╡ 2ed0a50b-2ab5-4bb7-9dbf-3af1f92c89c4
function chi_div_w(w::Vector, path::String, n::Int64)
	η = 0.001
	g = make_input(path) 
	gs = MasubaraGF(n, g.GF[length(g.GF)-n+1:end])
	return [spectrum_density(i,η,gs) for i in w]
end

# ╔═╡ 64b77123-bab8-48ed-9291-409f6158c8ba
begin
	n =  5
	ω = [ i for i in range(-10,10,step=1.e-2)]
		
	grange = [1.0, 1.5, 2.0]
	brange = [10, 20, 30, 40]
	chi =Matrix{Vector}(undef,4,3)
	for i = 1:4, j = 1:3
		path = @sprintf "../data/gt/gt28_g_%.1f_b_%i.txt" grange[j] brange[i]
		chi[i,j] = chi_div_w(ω, path, n)
	end
end

# ╔═╡ cc64bdcc-a924-43b6-afb3-adc491eb1d3a
begin
	plot(ω, chi[1,1], label="β = 10, g = 1.0")
	plot!(ω, chi[2,1], label="β = 20, g = 1.0")
	plot!(ω, chi[3,1], label="β = 30, g = 1.0")
	plot!(ω, chi[4,1], label="β = 40, g = 1.0")
	plot!(xlim=(-1,1))
end

# ╔═╡ 3fcbf32e-f5ab-4fc0-95a7-f40696ea5de6
csum = 0.01 .* sum.(chi)

# ╔═╡ 1aad6d9d-7943-4742-aead-e946f3cb0e42
# check sum rule

# ╔═╡ 6569707a-8c7c-4692-8162-07796a8c65f7
md"""
## Compare $G(\tau)$
"""

# ╔═╡ ac03e64f-4857-4634-a86d-bc351aa84d75
begin
	dtau_1p0_1 = readdlm("../data/TFIsing-Li/gtau_1.0_b20.txt")
	dtau_1p0_2 = readdlm("../data/gt/gf_t8_g_1.0_b_20.txt")
	dtau_1p0_3 = readdlm("../data/gt/gf_t28_g_1.0_b_20.txt")
	gtau_1p0 = zeros(401,3)
	gtau_1p0[:,1] = dtau_1p0_1[:,1]
	gtau_1p0[:,2] = dtau_1p0_2[:,2] .- dtau_1p0_1[:,2]
	gtau_1p0[:,3] = dtau_1p0_3[:,2] .- dtau_1p0_1[:,2]
end

# ╔═╡ 97ca7dd6-ea05-4f18-a9b2-d2ccca9f423f
begin
	plot(gtau_1p0[:,1], gtau_1p0[:,2], label="χ=8")
	plot!(gtau_1p0[:,1], gtau_1p0[:,3], label="χ=2 * 8")
	plot!(xlabel="τ/β", ylabel="G(τ) loss", title="g = 1.0, β=20")
end

# ╔═╡ 04674261-d673-4ccc-a3f6-fe87e6332375
begin
	dtau_1p5_1 = readdlm("../data/TFIsing-Li/gtau_1.5_b20.txt")
	dtau_1p5_2 = readdlm("../data/gt/gf_t8_g_1.5_b_20.txt")
	dtau_1p5_3 = readdlm("../data/gt/gf_t28_g_1.5_b_20.txt")
	gtau_1p5 = zeros(401,3)
	gtau_1p5[:,1] = dtau_1p5_1[:,1]
	gtau_1p5[:,2] = dtau_1p5_2[:,2] .- dtau_1p5_1[:,2]
	gtau_1p5[:,3] = dtau_1p5_3[:,2] .- dtau_1p5_1[:,2]
end

# ╔═╡ 976e673c-a58e-4f09-be4a-f09d3512ef23
begin
	plot(gtau_1p5[:,1], gtau_1p5[:,2], label="χ=8")
	plot!(gtau_1p5[:,1], gtau_1p5[:,3], label="χ=2 * 8")
	plot!(xlabel="τ/β", ylabel="G(τ) loss", title="g = 1.5, β=20")
end

# ╔═╡ eb93aae6-1b7d-4055-b35e-d3b543248a1a
begin
	dtau_2p0_1 = readdlm("../data/TFIsing-Li/gtau_2.0_b20.txt")
	dtau_2p0_2 = readdlm("../data/gt/gf_t8_g_2.0_b_20.txt")
	dtau_2p0_3 = readdlm("../data/gt/gf_t28_g_2.0_b_20.txt")
	gtau_2p0 = zeros(401,3)
	gtau_2p0[:,1] = dtau_2p0_1[:,1]
	gtau_2p0[:,2] = dtau_2p0_2[:,2] .- dtau_2p0_1[:,2]
	gtau_2p0[:,3] = dtau_2p0_3[:,2] .- dtau_2p0_1[:,2]
end

# ╔═╡ 515b2659-7624-44a7-8cfb-108fde0bf5fc
begin
	plot(gtau_2p0[:,1], gtau_2p0[:,2], label="χ=8")
	plot!(gtau_2p0[:,1], gtau_2p0[:,3], label="χ=2 * 8")
	plot!(xlabel="τ/β", ylabel="G(τ) loss", title="g = 2.0, β=20")
end

# ╔═╡ 524af3b4-05eb-4269-bd7e-3de6f2113783
begin
	g = 1.5
	a,b = divrem(g, 1.0)
	a = Integer(a)
	b = Integer(b * 10)
	
	"""
	#import data
	p = @sprintf "../data/TFIsing-Li/origional/g_%ip%i.txt" a b
	op = @sprintf "../data/TFIsing-Li/gt_%.1f_b20.txt" g
	
	d = readdlm(p) ; len = length(d)/2 |> Integer
	t = [i for i in range(0,200,step=0.5)]
	
	#export data [t re.gt im.gt]
	open(op,"w") do file
		for i = 1:len
			writedlm(file,[t[i] d[2*i-1] d[2*i]])
		end
	end
	"""
	
	"""
	p = @sprintf "../data/TFIsing-Li/origional/g_%ip%i_tau_N_501.txt" a b
	op = @sprintf "../data/TFIsing-Li/gtau_%.1f_b20.txt" g
	
	#import data
	d = readdlm(p) ; len = length(d)
	tau = [i for i in range(0,20,step=0.05)]
	
	#export data [τ/β val]
	open(op,"w") do file
		for i = 1:len
			writedlm(file,[tau[i]/20 d[i]])
		end
	end
	"""
	
	"""
	nval = "A"
	p = @sprintf "../data/TFIsing-Li/origional/g_%ip%i_%s.txt" a b nval
	op = @sprintf "../data/TFIsing-Li/%sw_%.1f_b20.txt" nval g
	
	#import data
	d = readdlm(p) ; len = length(d)
	
	#make omega
	o = [-i for i in range(-2π,2π,step = 4π/len)]
	
	#export data [ω val]
	open(op,"w") do file
		for i = 1:len
			writedlm(file,[o[i] d[i]])
		end
	end
	"""
	
	"reexport data"
end

# ╔═╡ 42720925-c2f4-4103-a0a3-1cb51bc680f5
function f(w::Real)
	res = 1 - exp(-20*w)
	return res
end

# ╔═╡ 4147ddc7-bc31-47aa-a02e-eb193382b83c
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

# ╔═╡ 0b403313-abc1-4e91-a1f2-c745f167740f
begin
	n1 =  20
	#Zi-Long Li data
	A1_li= readdlm("../data/TFIsing-Li/Aw_1.0_b20.txt")
	omega = A1_li[:,1]
	A1w_li = A1_li[:,2] ./ omega
	
	
	#theory sachdev
	A1_theory = [spectral_density(1.0, i, 20) for i in omega]
	
	#cmpo data
	g1 = make_input("../data/gt/gt28_g_1.0_b_20.txt") 
	g1s = MasubaraGF(10, g1.GF[length(g1.GF) - n1 + 1:end])
	A1_cmpo = [spectrum_density(w,η,g1s)*w for w in omega] 
	A1w_cmpo = [spectrum_density(w,η,g1s) for w in omega] 
end

# ╔═╡ 8060aa27-1de7-4b40-9665-15ce4760dee4
begin
	plot(omega, A1_theory, line=(2,:black), label="theory")
	plot!(omega, A1_li[:,2] ,line=(:dash,2), label="Li")
	plot!(omega, A1_cmpo, line=(:dash,2), label=@sprintf "χ=2*8,n=%i" n1)
	plot!(xlim=(0,5), ylim=(0,10))
	plot!(xlabel="ω", ylabel="A(ω)",title="g=1.0, β=20")
end

# ╔═╡ 1848b1c3-b52f-4b5a-9944-4d79e91a59a3
begin
	plot(omega, A1w_li,line=(2), label="Li")
	plot!(omega, A1w_cmpo, line=(:dash,2), label=@sprintf "χ=2*8,n=%i" n1)
	plot!(xlim=(-1,1))
	plot!(xlabel="ω", ylabel="A(ω)/ω",title="g=1.0, β=20")
end

# ╔═╡ 83ad4c7b-8174-4d63-bb3c-33aa4faa3edb
begin
	n1p5 =  5
	#Zi-Long Li data
	A1p5_li= readdlm("../data/TFIsing-Li/Aw_1.5_b20.txt")
	A1p5w_li = A1p5_li[:,2] ./ omega
	
	#cmpo data
	g1p5 = make_input("../data/gt/gt28_g_1.5_b_20.txt") 
	g1p5s = MasubaraGF(n1p5, g1p5.GF[length(g1p5.GF) - n1p5 + 1:end])
	A1p5_cmpo = [spectrum_density(w,η,g1p5s)*w for w in omega] 
	A1p5w_cmpo = [spectrum_density(w,η,g1p5s) for w in omega] 
end

# ╔═╡ fcb57525-938a-4951-8e21-11ca03a3b318
begin
	plot(omega, A1p5w_li ,line=(2), label="Li")
	plot!(omega, A1p5w_cmpo, line=(:dash,2), label=@sprintf "χ=2*8,n=%i" n1p5)
	plot!([-1,1], seriestype="vline",line=(:dash, :black))
	plot!(xlabel="ω", ylabel="A(ω)/ω",title="g=1.5, β=20")
end

# ╔═╡ 987bd147-f917-47a2-ad04-2571c9e69dd7
begin
	n2 =  5
	#Zi-Long Li data
	A2_li= readdlm("../data/TFIsing-Li/Aw_2.0_b20.txt")
	A2w_li = A2_li[:,2] ./ omega
	
	#cmpo data
	g2 = make_input("../data/gt/gt28_g_2.0_b_20.txt") 
	g2s = MasubaraGF(n2, g2.GF[length(g2.GF) - n2 + 1:end])
	A2_cmpo = [spectrum_density(w,η,g2s)*w for w in omega] 
	A2w_cmpo = [spectrum_density(w,η,g2s) for w in omega] 
end

# ╔═╡ cef09d6b-9fd8-47c8-9b6d-4def6c6d0a2a
begin
	plot(omega, A2w_li ,line=(2), label="Li")
	plot!(omega, A2w_cmpo, line=(:dash,2), label=@sprintf "χ=2*8,n=%i" n2)
	plot!([-2,2], seriestype="vline",line=(:dash, :black))
	plot!(xlabel="ω", ylabel="A(ω)/ω",title="g=1.5, β=20")
end

# ╔═╡ 73760863-b866-4f56-8115-8b03ab7c8192
pwd()

# ╔═╡ Cell order:
# ╟─93711dc9-3b99-4ce3-8b39-570b4a833f58
# ╟─25c4b6ca-c100-488e-8f73-e514f001aac8
# ╟─1ed9dc6a-2971-4c61-aa3c-b3f91f7cb34b
# ╟─0b403313-abc1-4e91-a1f2-c745f167740f
# ╟─8060aa27-1de7-4b40-9665-15ce4760dee4
# ╟─1848b1c3-b52f-4b5a-9944-4d79e91a59a3
# ╟─eaa15cb3-8958-4967-b3b1-51f509b83b8a
# ╟─83ad4c7b-8174-4d63-bb3c-33aa4faa3edb
# ╠═fcb57525-938a-4951-8e21-11ca03a3b318
# ╟─8f5723b4-5871-46b9-a706-495b1f34de77
# ╟─987bd147-f917-47a2-ad04-2571c9e69dd7
# ╟─cef09d6b-9fd8-47c8-9b6d-4def6c6d0a2a
# ╟─424f1e78-bbb1-4d55-8c14-a64d2f61ab0f
# ╟─2ed0a50b-2ab5-4bb7-9dbf-3af1f92c89c4
# ╟─64b77123-bab8-48ed-9291-409f6158c8ba
# ╟─cc64bdcc-a924-43b6-afb3-adc491eb1d3a
# ╠═3fcbf32e-f5ab-4fc0-95a7-f40696ea5de6
# ╠═1aad6d9d-7943-4742-aead-e946f3cb0e42
# ╟─6569707a-8c7c-4692-8162-07796a8c65f7
# ╟─ac03e64f-4857-4634-a86d-bc351aa84d75
# ╟─97ca7dd6-ea05-4f18-a9b2-d2ccca9f423f
# ╟─04674261-d673-4ccc-a3f6-fe87e6332375
# ╟─976e673c-a58e-4f09-be4a-f09d3512ef23
# ╟─eb93aae6-1b7d-4055-b35e-d3b543248a1a
# ╟─515b2659-7624-44a7-8cfb-108fde0bf5fc
# ╟─524af3b4-05eb-4269-bd7e-3de6f2113783
# ╟─42720925-c2f4-4103-a0a3-1cb51bc680f5
# ╟─4147ddc7-bc31-47aa-a02e-eb193382b83c
# ╟─73760863-b866-4f56-8115-8b03ab7c8192
# ╠═6d565956-314f-11ec-0005-a12d4cf1fb59
