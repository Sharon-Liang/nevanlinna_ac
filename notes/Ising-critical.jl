### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ fa4062f3-dd6c-4814-8b65-c6641b352345
begin
	using Pkg
	Pkg.activate("../")
	using nevanlinna_ac
	using Plots
	using DelimitedFiles
	using Printf
	"packages"
end

# ╔═╡ 7b97d05b-f414-47d4-b154-6a64a931c712
begin
	a = readdlm("../data/a_j_1.0.txt")
	#g =make_input("../data/g_j_1.0.txt")
	gn8 = make_input("../data/gn_1.0_chi8.txt")
	gn28 = make_input("../data/gn_1.0_chi28.txt")
	gt8 = make_input("../data/gt_1.0_chi8.txt")
	gt28 = make_input("../data/gt_1.0_chi28.txt")
	"a1 = A(ω), gn= Gn, gt = ∑_mn Cmn(iωn)/(Em - En)"
end

# ╔═╡ 9f67e944-d615-4346-b891-b0c73fc03289
function lessnode(n::Int, g::MasubaraGF)
	if n > g.number @error "n is too large" end
	b = g.number - n +1
	return MasubaraGF(n, g.GF[b:end])
end

# ╔═╡ bdf32191-0d89-45a2-9d96-f5886437f7d8
n=5

# ╔═╡ 9e672600-e6fa-4af7-b1e8-d955b417805f
begin
	ωn8 = [gt8.GF[i].ωn for i = 1:gt8.number]
	Gt8 = [gt8.GF[i].val |> imag for i = 1:gt8.number]
end

# ╔═╡ 8de07358-9370-4f89-a629-b2b5a75579f8
begin
	ωn28 = [gt28.GF[i].ωn for i = 1:gt28.number]
	Gt28 = [gt28.GF[i].val |> imag for i = 1:gt28.number]
end

# ╔═╡ 6b97e468-61cf-4ffe-8593-3169e67fa00f
begin
	scatter(ωn8, Gt8 .- Gt28, label="χ=8")
	plot!(xlabel="ωn",ylabel="∑_mn Cmn(iωn)/(Em - En)")
	#scatter!(ωn28, Gt28, label="χ=2*8")
end

# ╔═╡ 416469f4-38bf-4ae6-b3e0-07db08128cc8
omega = [i for i in range(0,5,length=1000)]

# ╔═╡ 276d16ed-e05c-455c-9b34-e54c383ecc19
η = 0.001

# ╔═╡ b5b2c40c-6aa2-4c03-87d0-87a3d91208cc
a1 = [spectrum_density(w,η,gt8)*w for w in omega] 

# ╔═╡ eedc42b3-f25a-4eca-bde6-f3c07b6c5bfd
a1b = [spectrum_density(w,η,gt8,etype=BigFloat)*w for w in omega] 

# ╔═╡ 84f9166b-6810-4828-823e-998bff6afee2
a2 = [spectrum_density(w,η,lessnode(n,gt8))*w for w in omega] 

# ╔═╡ 4ff996d4-f6cd-4860-861d-a44b6f39a502
a2b = [spectrum_density(w,η,lessnode(n,gt8),etype=BigFloat)*w for w in omega] 

# ╔═╡ f75579f2-09da-4999-9ca3-c97f720f3fff
a3 = [spectrum_density(w,η,gt28)*w for w in omega] 

# ╔═╡ 03ee5ceb-3a7e-4528-8433-90b5bf24d6a5
a3b = [spectrum_density(w,η,gt28,etype=BigFloat)*w for w in omega] 

# ╔═╡ 89919e17-9c3e-49bb-a64a-aed4443498d6
a4 = [spectrum_density(w,η,lessnode(n,gt28))*w for w in omega] 

# ╔═╡ 9042b672-4c9d-4b6a-8933-f25cb08944d1
begin
	plot(a[:,1], a[:,2], line=(:black,2), label="theory")
	#plot!(omega, a1,line=(:dash,2), label="N=30")
	#plot!(omega, a1b,line=(:dash,2), label="N=30, BigFloat")
	plot!(omega, a2,line=(2), label= @sprintf "N=%i" n)
	#plot!(omega, a2b,line=(:dash,2), label= @sprintf "N=%i, BigFloat" n)
	
	#plot!(omega, a3,line=(:dash,2), label="N=30,2*8")
	#plot!(omega, a3b,line=(:dash,2), label="N=30, 2*8, BigFloat")
	
	plot!(omega, a4,line=(2), label= @sprintf "N=%i, 2*8" n)
	#plot!(omega, a4b,line=(:dash,2), label= @sprintf "N=%i, 2*8, BigFloat" n)
	plot!(title="g=1.0, β=20")
	#plot!(ylim=(-0.1,10))
	plot!(xlabel="ω", ylabel="A(ω)",
		#xtickfont=font(12), 
    	#ytickfont=font(12), 
    	#guidefont=font(12), 
    	#legendfont=font(12)
		)
end

# ╔═╡ aa53e759-c0e4-4200-9031-3128b25e7f78
a4b = [spectrum_density(w,η,lessnode(n,gt28),etype=BigFloat)*w for w in omega] 

# ╔═╡ 1091ad4a-2a3a-11ec-0191-6bcfedbfd580
pwd()

# ╔═╡ Cell order:
# ╟─7b97d05b-f414-47d4-b154-6a64a931c712
# ╟─9f67e944-d615-4346-b891-b0c73fc03289
# ╟─bdf32191-0d89-45a2-9d96-f5886437f7d8
# ╟─b5b2c40c-6aa2-4c03-87d0-87a3d91208cc
# ╟─eedc42b3-f25a-4eca-bde6-f3c07b6c5bfd
# ╟─84f9166b-6810-4828-823e-998bff6afee2
# ╟─4ff996d4-f6cd-4860-861d-a44b6f39a502
# ╟─f75579f2-09da-4999-9ca3-c97f720f3fff
# ╟─03ee5ceb-3a7e-4528-8433-90b5bf24d6a5
# ╟─89919e17-9c3e-49bb-a64a-aed4443498d6
# ╟─aa53e759-c0e4-4200-9031-3128b25e7f78
# ╟─9042b672-4c9d-4b6a-8933-f25cb08944d1
# ╠═9e672600-e6fa-4af7-b1e8-d955b417805f
# ╟─8de07358-9370-4f89-a629-b2b5a75579f8
# ╟─6b97e468-61cf-4ffe-8593-3169e67fa00f
# ╟─416469f4-38bf-4ae6-b3e0-07db08128cc8
# ╟─276d16ed-e05c-455c-9b34-e54c383ecc19
# ╟─fa4062f3-dd6c-4814-8b65-c6641b352345
# ╟─1091ad4a-2a3a-11ec-0191-6bcfedbfd580
