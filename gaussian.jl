### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ a03fbe02-8cb1-4cfc-92f1-0bc2163a4160
begin
	using Pkg; Pkg.activate("./")
	using nevanlinna_ac
	using Plots
	using DelimitedFiles
	"packages"
end

# ╔═╡ 08862677-7f44-4325-ba74-b42d40e307ab
begin
	A1 = readdlm("./data/gaussian/A_delta_eta_0.05.txt")
	d1 = readdlm("./data/gaussian/giwn_delta_eta_0.05.txt")
	x1 = d1[:,1]
	y1 = d1[:,2] .+ 1.0im*d1[:,3]
	x1 = x1 |> reverse
	y1 = y1 |> reverse
	"delta function spectrum"
end

# ╔═╡ 859f5f5a-02fa-4efd-8f23-5d4a1860d008
num = 21

# ╔═╡ e491e00e-bbf3-47ee-92fa-b13518efe869
begin
	x1n, y1n = toNevanlinnadata(x1[1:num], y1[1:num],:f)
	isNevanlinnasolvable(x1n,y1n)
end

# ╔═╡ ee2bd9df-7355-4bce-b3a7-7a91efc4e91a
begin
	A1n, name = spectrum(A1[:,1],0.05,x1[1:num],y1[1:num],:f)
	plot(A1[:,1], A1[:,2], lw = 2, label="theory")
	plot!(A1[:,1], A1n, line=(:dash,2), label="nevanlinna")
	plot!(xlabel="ω", ylabel=name)
end

# ╔═╡ 87cc2e29-fe91-45c0-b609-385bff1756c8
begin
	A2 = readdlm("./data/gaussian/A_gaussian_1.txt")
	d2 = readdlm("./data/gaussian/giwn_gaussian_1.txt")
	x2 = d2[:,1]
	y2 = d2[:,2] .+ 1.0im*d2[:,3]
	x2 = x2 |> reverse
	y2 = y2 |> reverse
	"single gaussian function spectrum"
end

# ╔═╡ 02012e17-a7d3-47ef-b706-5bcd85978757
ng1 = 40

# ╔═╡ bca8dbc4-ab81-4f5c-9969-cde1c5ee981a
begin
	x2n, y2n = toNevanlinnadata(x2[1:ng1], y1[1:ng1],:f)
	isNevanlinnasolvable(x2n,y2n)
end

# ╔═╡ dea33c6f-8c10-4adc-928d-067ed97d3159
begin
	A2n, _ = spectrum(A2[:,1],0.05,x2[1:ng1],y2[1:ng1],:f)
	plot(A2[:,1], A2[:,2], lw = 2, label="theory")
	plot!(A2[:,1], A2n, line=(:dash,2), label="nevanlinna")
	plot!(xlabel="ω", ylabel=name)
end

# ╔═╡ 49517ab2-fb41-403d-8dde-0f47c56676e2
begin
	A3 = readdlm("./data/gaussian/A_gaussian_3.txt")
	d3 = readdlm("./data/gaussian/giwn_gaussian_3.txt")
	x3 = d3[:,1]
	y3 = d3[:,2] .+ 1.0im*d3[:,3]
	x3 = x3 |> reverse
	y3 = y3 |> reverse
	"multi gaussian function spectrum"
end

# ╔═╡ a8fa541d-d41d-4466-b8ae-2b4ca91b7193
ng2 = 40

# ╔═╡ 83631e80-e07b-483e-8b47-08928bee91af
begin
	x3n, y3n = toNevanlinnadata(x3[1:ng2], y3[1:ng2],:f)
	isNevanlinnasolvable(x3n,y3n)
end

# ╔═╡ 580e257c-9ccf-405f-b9ee-0db0470ebaa2
begin
	A3n, _ = spectrum(A3[:,1],0.05,x3[1:ng2],y3[1:ng2],:f)
	plot(A3[:,1], A3[:,2], lw = 2, label="theory")
	plot!(A3[:,1], A3n, line=(:dash,2), label="nevanlinna")
	plot!(xlabel="ω", ylabel=name)
end

# ╔═╡ 7b4d617c-62f9-11ec-2ea6-9937cf302a50
pwd()

# ╔═╡ Cell order:
# ╟─08862677-7f44-4325-ba74-b42d40e307ab
# ╟─859f5f5a-02fa-4efd-8f23-5d4a1860d008
# ╠═e491e00e-bbf3-47ee-92fa-b13518efe869
# ╠═ee2bd9df-7355-4bce-b3a7-7a91efc4e91a
# ╠═87cc2e29-fe91-45c0-b609-385bff1756c8
# ╠═02012e17-a7d3-47ef-b706-5bcd85978757
# ╠═bca8dbc4-ab81-4f5c-9969-cde1c5ee981a
# ╠═dea33c6f-8c10-4adc-928d-067ed97d3159
# ╟─49517ab2-fb41-403d-8dde-0f47c56676e2
# ╟─a8fa541d-d41d-4466-b8ae-2b4ca91b7193
# ╟─83631e80-e07b-483e-8b47-08928bee91af
# ╟─580e257c-9ccf-405f-b9ee-0db0470ebaa2
# ╟─a03fbe02-8cb1-4cfc-92f1-0bc2163a4160
# ╟─7b4d617c-62f9-11ec-2ea6-9937cf302a50
