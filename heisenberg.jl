### A Pluto.jl notebook ###
# v0.17.2

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

# ╔═╡ ed1194b4-5326-11ec-0b3a-e193af55a577
begin
	using Pkg
	Pkg.activate("./")
	using LinearAlgebra
	using Random; Random.seed!()
	using DelimitedFiles
	using Plots
	using StatsFuns, SpecialFunctions
	using nevanlinna_ac
	using Printf
	"packages"
end

# ╔═╡ 88a0ea81-b1d2-4d8f-beac-8cc0810f0111
omega = [i for i in range(-5,5,length=1600)]

# ╔═╡ 2d56c971-0e42-4e33-8d6d-2ca4ee612dd2
md"""
temperature range: $\beta = [1,2]$
"""

# ╔═╡ a15fbbb5-783a-4548-aeb6-39858fb9b17b
#@bind β html"<input type=range min=1 max=2 step=1>"

# ╔═╡ 3c9a2ed2-f839-4261-a151-3c0f5b1250bd
β =1

# ╔═╡ c2fc673e-edd2-47ba-bfe2-8a07e155b495
md"""
frequency number$= 1:40$
"""

# ╔═╡ 0fa0621d-c831-43f1-900c-c3cb8f4d4514
@bind freq_num html"<input type=range min=1 max=40 step=1>"

# ╔═╡ e8a66838-2d29-45e4-9097-cabe3ef689e8
begin
	pname1 = @sprintf "./data/xxz-imagtime/giwn_heisenberg_D_2m8_beta_%i.txt" β
	dh = readdlm(pname1);
	@sprintf "heisenberg imagtime β = %i: dh" β
end

# ╔═╡ 68af1400-f5d7-4269-8a52-ecd7844f9beb
begin
	pising = @sprintf "./data/ising-imagtime/giwn_g_1.0_D_2m8_beta_%i.txt" 10
	di = readdlm(pising)
	"ising imagtime β = 10: di"
end

# ╔═╡ 2644c7f4-3972-47d8-836c-8e5e758f3618
begin
	plot(di[:,1].*10/2π, di[:,3], line=(2), marker=(:circle),label="ising real")
	plot!(dh[:,1].*β/2π, dh[:,3], line=(:dash, 2), marker=(:circle),label="heisenberg real")
	plot!(xlabel="n", ylabel="real G(iωn)/ΔE",legend=:bottomright)
end

# ╔═╡ 5f494a87-13ac-4918-979f-5810b5eef03e
begin
	plot(di[:,1].*10, di[:,2], line=(2), marker=(:circle),label="ising imag")
	plot!(dh[:,1].*β, dh[:,2], line=(:dash, 2), marker=(:circle),label="heisenberg imag")
	plot!(xlabel="β*ωn", ylabel="imag G(iωn)",legend=:bottomright)
end

# ╔═╡ 75399062-1c82-48b3-82db-4f7bec7a7620
p1 = @sprintf "./data/xxz-imagtime/gtau_heisenberg_D_2m8_beta_%i.txt" β
p2 = @sprintf "./data/ising-imagtime/gtau_heisenberg_D_2m8_beta_%i.txt" 10
d1 = readdlm(p1); d2=readdlm(p2)
plot(d1[:,1]./10, d1[:,2], )

# ╔═╡ 87041b8b-d324-408a-bcc2-3c2cb640fdab
function chi_div_w(w::Vector, path::String, n::Int64)
	η = 0.001
	g = make_input(path) 
	gs = MasubaraGF(n, g.GF[length(g.GF)-n+1:end])
	return [spectrum_density(i,η,gs) for i in w]
end

# ╔═╡ 70c580d6-9730-4866-8676-6cdf1b7af54c
begin
	pname = @sprintf "./data/xxz-imagtime/gdivwn_heisenberg_D_2m8_beta_%i.txt" β
	AdivW = chi_div_w(omega, pname, freq_num)
	line_name = @sprintf "freq_number=%i" freq_num
	plot(omega, AdivW, lw=2, label=line_name)
	plot!(xlim=(-5,5), ylim=(-10,10))
	fig_name = @sprintf "AFM Heisenberg model: β=%i" β
	plot!(xlabel="ω", ylabel="A(ω)/ω", title=fig_name)
end

# ╔═╡ d12d9aed-7434-47bd-af5f-fc35ac059f89
begin
	plot(omega, AdivW .* omega ./(1 .- exp.(-β.*omega)), lw=2, label=line_name)
	plot!(xlim=(-5,5), ylim=(-10,10))
	plot!(xlabel="ω", ylabel="S(ω)", title=fig_name)
end

# ╔═╡ c9898fcb-5580-4214-be9b-2f1a5c18e856
function Masubara_freq(n::Int64, β::Real; type::Symbol=:b)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end

# ╔═╡ 284645b9-4e38-41c3-9592-5e9e0cc192c6
function giwn(n::Int, β::Real, τ::Vector, g::Vector)
    isapprox(β, τ[end], atol=0.1) ? iωn = 1.0im*Masubara_freq(n,β,type=:b) : error("wrong β")
    len = length(τ)
    res = 0.0 + 0.0im
    for i = 1: len
        res += -g[i] * exp(iωn*τ[i])
    end
    return β/len*res
end

# ╔═╡ Cell order:
# ╟─88a0ea81-b1d2-4d8f-beac-8cc0810f0111
# ╟─2d56c971-0e42-4e33-8d6d-2ca4ee612dd2
# ╠═a15fbbb5-783a-4548-aeb6-39858fb9b17b
# ╠═3c9a2ed2-f839-4261-a151-3c0f5b1250bd
# ╟─c2fc673e-edd2-47ba-bfe2-8a07e155b495
# ╟─0fa0621d-c831-43f1-900c-c3cb8f4d4514
# ╠═d12d9aed-7434-47bd-af5f-fc35ac059f89
# ╟─70c580d6-9730-4866-8676-6cdf1b7af54c
# ╠═e8a66838-2d29-45e4-9097-cabe3ef689e8
# ╠═68af1400-f5d7-4269-8a52-ecd7844f9beb
# ╠═2644c7f4-3972-47d8-836c-8e5e758f3618
# ╠═5f494a87-13ac-4918-979f-5810b5eef03e
# ╠═75399062-1c82-48b3-82db-4f7bec7a7620
# ╠═87041b8b-d324-408a-bcc2-3c2cb640fdab
# ╟─284645b9-4e38-41c3-9592-5e9e0cc192c6
# ╟─c9898fcb-5580-4214-be9b-2f1a5c18e856
# ╟─ed1194b4-5326-11ec-0b3a-e193af55a577
