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
	using PlutoUI
	"packages"
end

# ╔═╡ 88a0ea81-b1d2-4d8f-beac-8cc0810f0111
omega = [i for i in range(-5,5,length=1600)]

# ╔═╡ 90229d63-5e6f-465b-9717-c4a6a3f3cb53
beta = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0]

# ╔═╡ 2c351b10-52dd-442a-9741-c7c85961aa1b
η = 0.001

# ╔═╡ 32ef830a-922d-46f1-a973-9480a48fd5f3
@bind J Select([0.0, 2.0])

# ╔═╡ 46582e6c-00c9-4640-b322-8c658d6fd685
@bind func Select(["pm", "mp", "pz"])

# ╔═╡ a15fbbb5-783a-4548-aeb6-39858fb9b17b
@bind β Slider(beta)

# ╔═╡ 0fa0621d-c831-43f1-900c-c3cb8f4d4514
@bind freq_num html"<input type=range min=1 max=40 step=1>"

# ╔═╡ e8a66838-2d29-45e4-9097-cabe3ef689e8
begin
	p1 = @sprintf "./data/xxz/imagtime/giwn/Jz_%.1f_%s_D_8m2_beta_%i.txt" J func β
	x1, y1 = readGF(p1, num = freq_num )
	@sprintf "heisenberg imagtime β = %i " β
end

# ╔═╡ aadd636b-71dd-4e91-a981-a232fe29df72
begin
	x1n, y1n = toNevanlinnadata(x1, y1, :b)
	isNevanlinnasolvable(x1n, y1n)
end

# ╔═╡ 26002478-73bf-42dc-a011-0bbedb92b80b
begin
	A1, name = spectrum(omega, η, x1, y1, :b)
	plot(omega, A1, lw=2, label=@sprintf "n=%i" freq_num)
	fig_title = @sprintf "β = %i" β
	plot!(xlabel="ω", ylabel = name, title = fig_title)
end

# ╔═╡ c9898fcb-5580-4214-be9b-2f1a5c18e856
function Masubara_freq(n::Int64, β::Real; type::Symbol=:b)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end

# ╔═╡ Cell order:
# ╟─88a0ea81-b1d2-4d8f-beac-8cc0810f0111
# ╟─90229d63-5e6f-465b-9717-c4a6a3f3cb53
# ╟─2c351b10-52dd-442a-9741-c7c85961aa1b
# ╟─32ef830a-922d-46f1-a973-9480a48fd5f3
# ╠═46582e6c-00c9-4640-b322-8c658d6fd685
# ╟─a15fbbb5-783a-4548-aeb6-39858fb9b17b
# ╟─e8a66838-2d29-45e4-9097-cabe3ef689e8
# ╟─aadd636b-71dd-4e91-a981-a232fe29df72
# ╟─0fa0621d-c831-43f1-900c-c3cb8f4d4514
# ╟─26002478-73bf-42dc-a011-0bbedb92b80b
# ╟─c9898fcb-5580-4214-be9b-2f1a5c18e856
# ╟─ed1194b4-5326-11ec-0b3a-e193af55a577
