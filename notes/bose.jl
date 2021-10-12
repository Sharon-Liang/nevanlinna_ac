### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 1ad3919a-1f79-11ec-1a44-d55ae420c8a6
begin
	using Pkg
	Pkg.activate("../")
	using nevanlinna_ac
	using Plots
	using DelimitedFiles
	using Printf
	"packages"
end

# ╔═╡ 8f75c874-0834-48a2-93a7-5521177826e4
η = 0.001

# ╔═╡ 8153e87c-fa69-44e0-8895-11ef2ffddf7c
omega = [i for i in range(0,5,length=1000)]

# ╔═╡ ad7da214-6055-4d57-b010-7fb4dbe96ea6
md"""
### J = 0.0; Γ = 1.0
"""

# ╔═╡ 8027d1ea-1188-4fd3-b1b9-99a9640c59f4
begin
	a0 =readdlm("../data/a_j_0.0.txt")
	g0 =make_input("../data/g_j_0.0.txt")
	"a0, g0"
end

# ╔═╡ ab3829e7-4ffa-422b-a62c-6f99e08bcc1d
isNevanlinnasolvable(g0)

# ╔═╡ fdad034c-6a53-43d3-b888-e38e2291b328
begin
	A0 = [spectrum_density(w,η,g0) for w in omega] 
	g0m =make_input("../data/gm_j_0.0.txt")
	g0n =make_input("../data/gn_j_0.0.txt")
	A0mn = [spectrum_density(w,η,g0m) - spectrum_density(w,η,g0n) for w in omega] 

	plot(a0[:,1], a0[:,2], lw=2, label="theory")
	plot!(omega, A0,line=(:dash,2), label="interpolate")
	plot!(omega, A0mn,line=(:dash,2), label="interpolate-2")
end

# ╔═╡ 0f5abb05-38e7-4009-93b2-f2ff7a4669f8
begin
	plot([g0.GF[i].ωn for i = 1:g0.number],[g0.GF[i].val |> imag for i = 1:g0.number],
		line=(2), label="imag")
	plot!([g0.GF[i].ωn for i = 1:g0.number],[g0.GF[i].val |> real for i = 1:g0.number],
		line=(2), label="real")
	plot!(xlabel="ωn", ylabel="GF value", legend=:bottomright)
end

# ╔═╡ 16aca10e-f308-488b-9b40-9dff09885d55
md"""
### J = Γ = 1.0
"""

# ╔═╡ bc61d277-7b0f-459f-bce9-8f0573a09f5b
begin
	a1 =readdlm("../data/a_j_1.0.txt")
	g1 =make_input("../data/g_j_1.0.txt")
	"a1 = A(ω), g1 = GF"
end

# ╔═╡ da4f1bcb-2ece-49ad-b9de-adf9e5d01e67
begin
	plot([g1.GF[i].ωn for i = 1:g1.number],[g1.GF[i].val |> imag for i = 1:g1.number],
		line=(2), label="imag")
	plot!([g1.GF[i].ωn for i = 1:g1.number],[g1.GF[i].val |> real for i = 1:g1.number],
		line=(2), label="real")
	plot!(xlabel="ωn", ylabel="GF value", legend=:bottomright)
end

# ╔═╡ 7b349ab2-f729-4a5c-99ae-6384ed3b3867
isNevanlinnasolvable(g1)

# ╔═╡ 1224a52a-55d0-4382-b678-544aa5ffe73d
A1 = [spectrum_density(w,η,g1) for w in omega]

# ╔═╡ 94d0a88f-f491-40c5-87f1-f4138f284fe3
begin
	g1m =make_input("../data/gm_j_1.0.txt")
	#g1n =make_input("../data/gn_j_1.0.txt")
	"∑_mn Cmn(iωn)/(Em - En)"
end

# ╔═╡ 43de4aed-4d19-4a8c-85f9-2ca5cd2cd2a8
begin
	plot([g1m.GF[i].ωn for i = 1:g1m.number],[g1m.GF[i].val |> imag for i = 1:g1m.number],
		line=(2), label="imag")
	plot!([g1m.GF[i].ωn for i = 1:g1m.number],[g1m.GF[i].val |> real for i = 1:g1m.number],
		line=(2), label="real")
	plot!(xlabel="ωn", ylabel="GF value", legend=:bottomright)
end

# ╔═╡ 96788e8e-0851-4a74-9214-9e6f7ae01852
A1m = [spectrum_density(w,η,g1m)*w for w in omega] 

# ╔═╡ 54064a72-1a7b-4c7d-99f6-5ce292df7970
plot(omega, A1m)

# ╔═╡ 9415fb29-4294-47ec-94ed-28d778c504aa
begin
	n = 5
	b = 30 - n + 1
	g1ms = MasubaraGF(n, g1m.GF[b:30])
	A1ms = [spectrum_density(w,η,g1ms)*w for w in omega] 
end

# ╔═╡ 1ab3f1f3-2c74-4717-91a7-f77863945d69
plot(omega, A1ms)

# ╔═╡ bd45648e-e070-40d2-b4e9-6b16f2c965af
isNevanlinnasolvable(g1m)

# ╔═╡ e150dbde-effb-417f-91ab-0ecc18d09945
function spectrum(ω::Real, η::Real, G::MasubaraGF)
    z = ω + 1.0im * η
    return -2im*evaluation(z, G) |> imag
end

# ╔═╡ 765f4fb6-43cb-44dc-b37a-e4eaf7e18197
begin
	ig1 =make_input("../data/ig_j_1.0.txt")
	"../data/ig_j_1.0.txt"
end

# ╔═╡ 75ed89a9-01f7-4587-bd21-c0ced4f8b819
ig1s = MasubaraGF(4, ig1.GF[27:30])

# ╔═╡ 13d87ae5-923c-476b-84b2-56bb44e46713
begin
	plot([ig1.GF[i].ωn for i = 1:ig1.number],[ig1.GF[i].val |> imag for i=1:ig1.number],
		line=(2), label="imag")
	plot!([ig1.GF[i].ωn for i = 1:ig1.number],
		  [ig1.GF[i].val |> real for i=1:ig1.number], 
		  line=(2), label="real")
	plot!(xlabel="ωn", ylabel="iGF value", legend=:bottomright)
end

# ╔═╡ a6049bf4-715f-4507-a43e-4566d1fa6dc9
isNevanlinnasolvable(ig1)

# ╔═╡ 59a7814e-d8ea-436d-ae13-5139c1dea323
isNevanlinnasolvable(ig1s)

# ╔═╡ 23de7386-30c1-406f-89bc-d9cb3ae47868
begin
	g1n =make_input("../data/gn_j_1.0.txt")
	"iωn * G(iωn)"
end

# ╔═╡ 85e0fbba-80e2-4f14-836a-3579e23becaa
isNevanlinnasolvable(g1n)

# ╔═╡ 8ec7f3b6-b45e-448e-b7c6-fe12c9309a2b
begin
	plot([g1n.GF[i].ωn for i = 1:g1n.number],[g1n.GF[i].val |> imag for i = 1:g1n.number],
		line=(2), label="imag")
	plot!([g1n.GF[i].ωn for i = 1:g1n.number],[g1n.GF[i].val |> real for i = 1:g1n.number],
		line=(2), label="real")
	plot!(xlabel="ωn", ylabel="GF value", legend=:bottomright)
end

# ╔═╡ d710e801-cb4b-4ff5-a125-13578c1fafb5
A1ns = [spectrum_density(w,η,g1n) for w in omega] 

# ╔═╡ 758d84d0-f844-4c8a-9c49-f6ae20eec7b6
begin
	
	#plot!(omega, A1,line=(:dash,2), label="interpolate")
	#plot!(omega, A1m,line=(2), label="C/dE")
	plot(omega, A1ms,line=(3), label="interpolate")
	plot!(a1[:,1], a1[:,2], line=(:black,:dash,3), label="theory")
	#plot!(omega, A1n,line=(:dash,2), label="C*ω")
	plot!(title="g=1.0, β=20")
	plot!(ylim=(-0.1,10))
	plot!(xlabel="ω", ylabel="A(ω)",
	xtickfont=font(12), 
    ytickfont=font(12), 
    guidefont=font(12), 
    legendfont=font(12))
end

# ╔═╡ 02fb2a24-5aa8-4e6e-9bff-68cefaa21f3b
md"""
$A(x) = \frac{x}{√(2π)} e^{-x^2 /2 }$
"""

# ╔═╡ 11f6a78a-4bf1-4a12-8176-43bcd9f1516a
begin
	A = readdlm("../data/ag1_bose.txt")
	"../data/ag1_bose.txt"
end

# ╔═╡ 791a5b9b-b2b3-4dcb-96c6-9d4ff5a0b083
begin
	g = make_input("../data/ig1_bose.txt")
	"../data/ig1_bose.txt"
end

# ╔═╡ 54a57c63-ff28-44b6-946c-6baeb5b5c56e
isNevanlinnasolvable(g1)

# ╔═╡ c243728f-64d2-4994-b686-44ef14b5dbaa
a = [spectrum_density(w,η,g)*w for w in omega]

# ╔═╡ 1600847a-de36-4c4a-bf74-d868b3bedaad
begin
	plot(A[:,1], A[:,2], lw=2, label="theory")
	plot!(omega, a,line=(:dash,2), label="interpolate")
	plot!(xlim=(0,5))
end

# ╔═╡ 286f864c-ab38-4ae0-8be6-66dea42c0056
begin
	plot([g.GF[i].ωn for i = 1:g.number],[g.GF[i].val |> imag for i = 1:g.number],
		line=(2), label="imag")
	plot!([g.GF[i].ωn for i = 1:g.number],[g.GF[i].val |> real for i = 1:g.number],
		line=(2), label="real")
	plot!(xlabel="ωn", ylabel="GF value", legend=:bottomright)
end

# ╔═╡ 84a0e459-36c9-462e-bd3c-14f24cb486c0
md"""
$A(x) = \frac{x^3}{√(2π)} e^{-x^2 /2 }$
"""

# ╔═╡ 3ca7de5c-9a31-41db-bad2-90b46c7111c5
begin
	A3 = readdlm("../data/ag3_bose.txt")
	"../data/ag3_bose.txt"
end

# ╔═╡ f21b14b4-fab8-40db-b76c-015f9f3c9e8d
begin
	g3 = make_input("../data/ig3_bose.txt")
	"../data/ig3_bose.txt"
end

# ╔═╡ f8e48a07-5d9f-486c-8599-08ea50747508
isNevanlinnasolvable(g3)

# ╔═╡ 9a3aeef6-d92f-49ee-ba9b-90e9f79d648b
a3 = [spectrum_density(w,η,g3)*w for w in omega]

# ╔═╡ 423c967d-3cfc-4650-9a1a-6d5dffedec0a
begin
	plot(A3[:,1], A3[:,2], lw=2, label="theory")
	plot!(omega, a3,line=(:dash,2), label="interpolate")
	plot!(xlim=(0,5))
end

# ╔═╡ d46b79d4-bfc7-4f85-a801-409e414e8439
begin
	plot([g3.GF[i].ωn for i = 1:g3.number],[g3.GF[i].val |> imag for i = 1:g3.number],
		line=(2), label="imag")
	plot!([g3.GF[i].ωn for i = 1:g3.number],[g3.GF[i].val |> real for i = 1:g3.number],
		line=(2), label="real")
	plot!(xlabel="ωn", ylabel="GF value", legend=:bottomright)
end

# ╔═╡ d44b60bb-6cff-4e81-9682-7bcf5643717e
pwd()

# ╔═╡ Cell order:
# ╟─8f75c874-0834-48a2-93a7-5521177826e4
# ╠═8153e87c-fa69-44e0-8895-11ef2ffddf7c
# ╟─ad7da214-6055-4d57-b010-7fb4dbe96ea6
# ╟─8027d1ea-1188-4fd3-b1b9-99a9640c59f4
# ╟─ab3829e7-4ffa-422b-a62c-6f99e08bcc1d
# ╠═fdad034c-6a53-43d3-b888-e38e2291b328
# ╟─0f5abb05-38e7-4009-93b2-f2ff7a4669f8
# ╟─16aca10e-f308-488b-9b40-9dff09885d55
# ╠═bc61d277-7b0f-459f-bce9-8f0573a09f5b
# ╟─da4f1bcb-2ece-49ad-b9de-adf9e5d01e67
# ╟─7b349ab2-f729-4a5c-99ae-6384ed3b3867
# ╠═1224a52a-55d0-4382-b678-544aa5ffe73d
# ╠═94d0a88f-f491-40c5-87f1-f4138f284fe3
# ╟─43de4aed-4d19-4a8c-85f9-2ca5cd2cd2a8
# ╠═96788e8e-0851-4a74-9214-9e6f7ae01852
# ╠═54064a72-1a7b-4c7d-99f6-5ce292df7970
# ╠═9415fb29-4294-47ec-94ed-28d778c504aa
# ╠═1ab3f1f3-2c74-4717-91a7-f77863945d69
# ╟─bd45648e-e070-40d2-b4e9-6b16f2c965af
# ╠═e150dbde-effb-417f-91ab-0ecc18d09945
# ╠═765f4fb6-43cb-44dc-b37a-e4eaf7e18197
# ╠═75ed89a9-01f7-4587-bd21-c0ced4f8b819
# ╟─13d87ae5-923c-476b-84b2-56bb44e46713
# ╠═85e0fbba-80e2-4f14-836a-3579e23becaa
# ╟─a6049bf4-715f-4507-a43e-4566d1fa6dc9
# ╟─59a7814e-d8ea-436d-ae13-5139c1dea323
# ╠═23de7386-30c1-406f-89bc-d9cb3ae47868
# ╟─8ec7f3b6-b45e-448e-b7c6-fe12c9309a2b
# ╠═d710e801-cb4b-4ff5-a125-13578c1fafb5
# ╠═758d84d0-f844-4c8a-9c49-f6ae20eec7b6
# ╟─02fb2a24-5aa8-4e6e-9bff-68cefaa21f3b
# ╟─11f6a78a-4bf1-4a12-8176-43bcd9f1516a
# ╟─791a5b9b-b2b3-4dcb-96c6-9d4ff5a0b083
# ╟─54a57c63-ff28-44b6-946c-6baeb5b5c56e
# ╠═c243728f-64d2-4994-b686-44ef14b5dbaa
# ╠═1600847a-de36-4c4a-bf74-d868b3bedaad
# ╠═286f864c-ab38-4ae0-8be6-66dea42c0056
# ╟─84a0e459-36c9-462e-bd3c-14f24cb486c0
# ╟─3ca7de5c-9a31-41db-bad2-90b46c7111c5
# ╟─f21b14b4-fab8-40db-b76c-015f9f3c9e8d
# ╟─f8e48a07-5d9f-486c-8599-08ea50747508
# ╟─9a3aeef6-d92f-49ee-ba9b-90e9f79d648b
# ╠═423c967d-3cfc-4650-9a1a-6d5dffedec0a
# ╟─d46b79d4-bfc7-4f85-a801-409e414e8439
# ╠═1ad3919a-1f79-11ec-1a44-d55ae420c8a6
# ╠═d44b60bb-6cff-4e81-9682-7bcf5643717e
