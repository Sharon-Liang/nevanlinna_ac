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

# ╔═╡ 3fa261a6-5334-11ec-1dd7-8b427445a7ab
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

# ╔═╡ c36eb877-d015-4947-8d95-61538a636646
g = 1.0

# ╔═╡ 6b130640-5f2b-49fb-8520-968f840af084
D = 8

# ╔═╡ ec19d6ec-2a4e-468b-8b87-d45a25d9a747
#beta = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0]
beta = [10.0, 20.0, 30.0, 40.0]

# ╔═╡ 0be8ccc2-00ab-46c6-8aa9-f7173d2c89cd
omega = [i for i in range(-4π,4π,step = π/400)]

# ╔═╡ 6f79f431-c46b-4ee0-be92-efaf0ace046f
η = 0.001

# ╔═╡ 8c2bd68e-1ba4-46e1-a07c-a749e1e4d1dc
@bind β Slider(beta)

# ╔═╡ 80060631-07fd-45f6-b0e2-585e65325c0c
begin
	pl = @sprintf "./data/ising/spectrum-Li/Aw/g_%.1f_beta_%i.txt" g β
	Al = readdlm(pl)
	"Zilong-Li data"
end

# ╔═╡ 02674a34-9dab-4431-bb18-6e3d7e7dc61b
begin
    p1 = @sprintf "./data/ising/imagtime/giwn/g_%.1f_D_%i_beta_%i.txt" g D β
    p2 = @sprintf "./data/ising/imagtime/giwn/g_%.1f_D_%im2_beta_%i.txt" g D β
	"G(iωn) data path"
end

# ╔═╡ 3fc73642-444a-47ef-9fd4-86ff786611df
begin
    p3 = @sprintf "./data/ising/imagtime/gpdE/g_%.1f_D_%i_beta_%i.txt" g D β
    p4 = @sprintf "./data/ising/imagtime/gpdE/g_%.1f_D_%im2_beta_%i.txt" g D β
	"G(iωn)* dE data path"
end

# ╔═╡ a44c4474-f9ac-4b6a-a861-19c786735642
@bind n1 html"<input type=range min=1 max=40 step=1>"

# ╔═╡ 8e9b19d0-ba7d-4ce4-8eca-1dd751826b56
begin
	x1, y1 = readGF(p1, num=n1)
	x1n, y1n = toNevanlinnadata(x1, y1, :b)
	isNevanlinnasolvable(x1n, y1n)
end

# ╔═╡ fe044405-f0f0-4425-8ea2-7320224ee62c
begin
	A1, name = spectrum(omega, η, x1, y1, :b)
	#A2, _ = spectrum(omega, η, x1, y1, :b)
	name
end

# ╔═╡ 9abcefdc-63b8-40ac-8fdd-d73957a7dc7e
begin
	x2, y2 = readGF(p2, num=n1)
	x2n, y2n = toNevanlinnadata(x2, y2, :b)
	isNevanlinnasolvable(x2n, y2n)
end

# ╔═╡ b7941790-0cbd-4c31-9643-4cffbc2c12ab
begin
	x3, y3 = readGF(p3, num = n1)
	x4, y4 = readGF(p4, num = n1)
	"G(iωn)* dE data"
end

# ╔═╡ fb53bede-302b-4161-a631-9b365db8140d
eltype(A1)

# ╔═╡ 12ddbf72-8c0a-4a95-9fb6-64c889ccf81b
begin
	x = copy(x1n)
	y = mt.(copy(y1n), 1.0im)
	"Data for generalized schure algorithm"
end

# ╔═╡ 23ac2281-4353-47e5-909c-f2884e589592
begin
	θ = [i for i in range(0,2π,length=500)]
	plot(cos.(θ), sin.(θ), line=(:black, 2),label=false)
	scatter!(real.(y), imag.(y), label="target data")
	plot!(title="target data")
end

# ╔═╡ cb20df42-476e-45f7-9fa6-e0e4b7b34e98
z = omega .+ 1.0im * η

# ╔═╡ da51c81d-5c82-454b-b999-3464c51e81d6
res = generalized_schur(z, x, y)

# ╔═╡ d8091d31-9ee6-4cea-b31b-e25bdfc6ee56
sparam = schur_parameter(x, y)

# ╔═╡ 9315b25c-c546-4a4e-bd96-89e4f72c6bf5
plot(omega, 1.0 .- abs.(res))

# ╔═╡ ded77134-54d7-4df9-b573-1830cb0b5295
@bind xlimit2 html"<input type=range min=0.1 max=5 step=1.e-2>"

# ╔═╡ 2013cccc-e046-4775-8344-4964258f049b
begin
	function mti(z::Number) 
	    return mt(z, 1.0im)
	end
	
	function imti(z::Number)
	    return imt(z, 1.0im)
	end
end

# ╔═╡ cb65501f-9af9-49ac-9aa1-ba76859aa202
begin
	scatter(imag(x1n), imag(y1n), label="imag -iωnG(iωn)")
	scatter!(x3, imag(y3), marker=(:star), label="imag G(iωn) * dE")
	scatter!(imag(x), imag.(imti.(generalized_schur(x, x, y))), marker=(:hexagon), label="-zG(z)")
end

# ╔═╡ f581e951-75ca-402a-876e-83c44939737c
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

# ╔═╡ d2e7f35f-5ef1-402c-b4ed-62682515edc9
begin
	de = [spectral_density(1.0, i, β) for i in omega]
	plot(omega, de .* omega, line=(:black,2),label="theory")
	plot!(Al[:,1], Al[:,2] .* Al[:,1], line=(:solid, 2),label="Li")
	label12 = @sprintf "χ=%i" D
	label22 = @sprintf "χ=%i×2" D
	plot!(omega, A1, line=(:dash, 2),label=label12)
	#plot!(omega, A2, line=(:dash, 2),label=label22)
	#plot!(xlim=(0, xlimit2), ylim=(0,maximum(de)*1.3))
	fig_title2=@sprintf "ising D=%i, β=%i, freq_num=%i" D β n1
	plot!(xlabel="ω", ylabel=name, title=fig_title2)
end

# ╔═╡ 7b3e2635-5a3f-4a65-b8bb-f2170155677b
function Masubara_freq(n::Int64, β::Real; type::Symbol=:b)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end

# ╔═╡ 0c6f07df-9192-4d31-bb03-271fcf56864f
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
# ╟─c36eb877-d015-4947-8d95-61538a636646
# ╟─6b130640-5f2b-49fb-8520-968f840af084
# ╟─ec19d6ec-2a4e-468b-8b87-d45a25d9a747
# ╟─0be8ccc2-00ab-46c6-8aa9-f7173d2c89cd
# ╟─6f79f431-c46b-4ee0-be92-efaf0ace046f
# ╟─80060631-07fd-45f6-b0e2-585e65325c0c
# ╟─02674a34-9dab-4431-bb18-6e3d7e7dc61b
# ╟─3fc73642-444a-47ef-9fd4-86ff786611df
# ╟─8e9b19d0-ba7d-4ce4-8eca-1dd751826b56
# ╟─9abcefdc-63b8-40ac-8fdd-d73957a7dc7e
# ╟─b7941790-0cbd-4c31-9643-4cffbc2c12ab
# ╟─cb65501f-9af9-49ac-9aa1-ba76859aa202
# ╟─fe044405-f0f0-4425-8ea2-7320224ee62c
# ╠═8c2bd68e-1ba4-46e1-a07c-a749e1e4d1dc
# ╠═a44c4474-f9ac-4b6a-a861-19c786735642
# ╟─d2e7f35f-5ef1-402c-b4ed-62682515edc9
# ╟─fb53bede-302b-4161-a631-9b365db8140d
# ╠═12ddbf72-8c0a-4a95-9fb6-64c889ccf81b
# ╠═23ac2281-4353-47e5-909c-f2884e589592
# ╠═cb20df42-476e-45f7-9fa6-e0e4b7b34e98
# ╠═da51c81d-5c82-454b-b999-3464c51e81d6
# ╟─d8091d31-9ee6-4cea-b31b-e25bdfc6ee56
# ╠═9315b25c-c546-4a4e-bd96-89e4f72c6bf5
# ╠═ded77134-54d7-4df9-b573-1830cb0b5295
# ╟─2013cccc-e046-4775-8344-4964258f049b
# ╟─f581e951-75ca-402a-876e-83c44939737c
# ╟─0c6f07df-9192-4d31-bb03-271fcf56864f
# ╟─7b3e2635-5a3f-4a65-b8bb-f2170155677b
# ╟─3fa261a6-5334-11ec-1dd7-8b427445a7ab
