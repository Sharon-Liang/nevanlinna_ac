### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a03fbe02-8cb1-4cfc-92f1-0bc2163a4160
begin
	using Pkg; Pkg.activate("./")
	using nevanlinna_ac
	using Plots
	using DelimitedFiles
	using Printf
	"packages"
end

# ╔═╡ 4912866a-e93f-4c52-971f-677d4dab2a26
@bind n1 html"<input type=range min=1 max=40 step=1>"

# ╔═╡ 08862677-7f44-4325-ba74-b42d40e307ab
begin
	A1 = readdlm("./data/gaussian/A_delta_eta_0.05.txt")
	x1, y1 = readGF("./data/gaussian/giwn_delta_eta_0.05.txt", rev=false, num = n1)
	"delta function spectrum"
end

# ╔═╡ e491e00e-bbf3-47ee-92fa-b13518efe869
begin
	x1n, y1n = toNevanlinnadata(x1, y1, :f)
	isNevanlinnasolvable(x1n,y1n)
end

# ╔═╡ ee2bd9df-7355-4bce-b3a7-7a91efc4e91a
begin
	A1n, name = spectrum(A1[:,1],0.05,x1,y1,:f)
	plot(A1[:,1], A1[:,2], lw = 2, label="theory")
	plot!(A1[:,1], A1n, line=(:dash,2), label=@sprintf "nevanlinna, n=%i" n1)
	plot!(xlabel="ω", ylabel=name)
end

# ╔═╡ f10b4706-8c54-4549-a38a-63a2c8ef22d3
@bind n2 html"<input type=range min=1 max=40 step=1>"

# ╔═╡ 87cc2e29-fe91-45c0-b609-385bff1756c8
begin
	A2 = readdlm("./data/gaussian/A_gaussian_1.txt")
	x2, y2 = readGF("./data/gaussian/giwn_gaussian_1.txt", num=n2)
	"single gaussian function spectrum"
end

# ╔═╡ bca8dbc4-ab81-4f5c-9969-cde1c5ee981a
begin
	x2n, y2n = toNevanlinnadata(x2, y2 ,:f)
	isNevanlinnasolvable(x2n,y2n)
end

# ╔═╡ dea33c6f-8c10-4adc-928d-067ed97d3159
begin
	A2n, _ = spectrum(A2[:,1],0.05,x2,y2,:f)
	plot(A2[:,1], A2[:,2], lw = 2, label="theory")
	plot!(A2[:,1], A2n, line=(:dash,2), label=@sprintf "nevanlinna, n=%i" n2)
	plot!(xlabel="ω", ylabel=name)
end

# ╔═╡ aaa9e66c-40f1-4f66-bf1f-2bf4b134cfe0
@bind n3 html"<input type=range min=1 max=40 step=1>"

# ╔═╡ 49517ab2-fb41-403d-8dde-0f47c56676e2
begin
	A3 = readdlm("./data/gaussian/A_gaussian_3.txt")
	x3, y3 = readGF("./data/gaussian/giwn_gaussian_3.txt", rev=true, num=n3)
	"multi gaussian function spectrum"
end

# ╔═╡ 83631e80-e07b-483e-8b47-08928bee91af
begin
	x3n, y3n = toNevanlinnadata(x3, y3,:f)
	isNevanlinnasolvable(x3n,y3n)
end

# ╔═╡ 580e257c-9ccf-405f-b9ee-0db0470ebaa2
begin
	A3n, _ = spectrum(A3[:,1],0.05,x3,y3,:f)
	plot(A3[:,1], A3[:,2], lw = 2, label="theory")
	plot!(A3[:,1], A3n, line=(:dash,2), label=@sprintf "nevanlinna, n=%i" n3)
	plot!(xlabel="ω", ylabel=name)
end

# ╔═╡ 7b4d617c-62f9-11ec-2ea6-9937cf302a50
pwd()

# ╔═╡ Cell order:
# ╠═08862677-7f44-4325-ba74-b42d40e307ab
# ╟─e491e00e-bbf3-47ee-92fa-b13518efe869
# ╠═4912866a-e93f-4c52-971f-677d4dab2a26
# ╠═ee2bd9df-7355-4bce-b3a7-7a91efc4e91a
# ╟─87cc2e29-fe91-45c0-b609-385bff1756c8
# ╟─bca8dbc4-ab81-4f5c-9969-cde1c5ee981a
# ╠═f10b4706-8c54-4549-a38a-63a2c8ef22d3
# ╠═dea33c6f-8c10-4adc-928d-067ed97d3159
# ╠═49517ab2-fb41-403d-8dde-0f47c56676e2
# ╟─83631e80-e07b-483e-8b47-08928bee91af
# ╟─aaa9e66c-40f1-4f66-bf1f-2bf4b134cfe0
# ╟─580e257c-9ccf-405f-b9ee-0db0470ebaa2
# ╟─a03fbe02-8cb1-4cfc-92f1-0bc2163a4160
# ╟─7b4d617c-62f9-11ec-2ea6-9937cf302a50
