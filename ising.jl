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
	"packages"
end

# ╔═╡ 657a6556-a153-4a62-8294-a157401cfca5
md"""
bond dimension: $D = [8, 16]$
"""

# ╔═╡ 6b130640-5f2b-49fb-8520-968f840af084
@bind D html"<input type=range min=8 max=16 step=8>"

# ╔═╡ 547b3223-0cf5-4331-a5ef-0fe9b64bd16d
md"""
temperature range: $\beta = [10,20,30,40]$
"""

# ╔═╡ 8c2bd68e-1ba4-46e1-a07c-a749e1e4d1dc
@bind β html"<input type=range min=10 max=40 step=10>"

# ╔═╡ 75fdb0f0-df02-4567-80f0-165c5e552d99
md"""
frequency number$= 1:40$
"""

# ╔═╡ a44c4474-f9ac-4b6a-a861-19c786735642
@bind freq_num html"<input type=range min=1 max=40 step=1>"

# ╔═╡ c7972842-afb9-4440-b35f-a7ffe45fa5b7
@bind xlimit html"<input type=range min=0.1 max=3 step=1.e-2>"

# ╔═╡ ded77134-54d7-4df9-b573-1830cb0b5295
@bind xlimit2 html"<input type=range min=0.1 max=5 step=1.e-2>"

# ╔═╡ c36eb877-d015-4947-8d95-61538a636646
g = 1.0

# ╔═╡ 7a9d8aa3-fd9d-43d0-bec6-92dcb59df920
(d,r) = divrem(g, 1);

# ╔═╡ ec19d6ec-2a4e-468b-8b87-d45a25d9a747
invT = [10, 20, 30, 40]

# ╔═╡ 80060631-07fd-45f6-b0e2-585e65325c0c
begin
	plca = Vector{String}(undef, length(invT))
	plcs = Vector{String}(undef, length(invT))
	for i = 1:length(invT)
	    plca[i] = @sprintf "./data/ising-Li/g_%ip%i_beta_%i_A.txt" d 10*r invT[i]
	    plcs[i] = @sprintf "./data/ising-Li/g_%ip%i_beta_%i_S.txt" d 10*r invT[i]
	end
	"Zilong-Li data"
end

# ╔═╡ 0be8ccc2-00ab-46c6-8aa9-f7173d2c89cd
omega = [i for i in range(-4π,4π,step = π/400)]

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

# ╔═╡ 8881d4a0-8bdf-4538-89a5-110d92c4f759
function chi_div_w(path::String)
    dl = readdlm(path)
    dl[:,2] = dl[:,2] ./ dl[:,1]
    return dl
end

# ╔═╡ 3d7f5cf2-6a6e-4900-b4f4-bdfab9997dcf
function chi_div_w(w::Vector, path::String, n::Int64)
	η = 0.001
	g = make_input(path) 
	gs = MasubaraGF(n, g.GF[length(g.GF)-n+1:end])
	return [spectrum_density(i,η,gs) for i in w]
end

# ╔═╡ 02674a34-9dab-4431-bb18-6e3d7e7dc61b
begin
	pc = Matrix{String}(undef, length(invT), 2)
	for i = 1:length(invT)
    	pc[i,1] = @sprintf "./data/ising-imagtime/gdivwn_g_%.1f_D_%i_beta_%i.txt" g D invT[i]
    	pc[i,2] = @sprintf "./data/ising-imagtime/gdivwn_g_%.1f_D_2m%i_beta_%i.txt" g D invT[i]
	end
	dc = Matrix{Vector}(undef, length(invT), 2)
	dl = Vector{Matrix}(undef,length(invT))
	for i = 1:length(invT)
    	dl[i] = chi_div_w(plca[i])
    	d = readdlm(plcs[i]); bs0 = invT[i]*d[800,2]
    	dl[i][800,2] = bs0
    	for j=1:2
        	dc[i,j] = chi_div_w(omega, pc[i,j], freq_num)
    	end
	end
end

# ╔═╡ 2582d5b5-1bbd-4e68-8324-77f26af0116d
begin
	t = div(β,10)
	de = [spectral_density(1.0, i, β) for i in omega]
	plot(omega, de ./omega, line=(:black,2),label="theory")
	plot!(dl[t][:,1], dl[t][:,2], line=(:solid, 2),label="Li")
	label1 = @sprintf "χ=%i" D
	label2 = @sprintf "χ=%i×2" D
	plot!(omega, dc[t,1],line=(:dash, 2),label=label1)
	plot!(omega, dc[t,2],line=(:dash, 2),label=label2)
	plot!(xlim=(-xlimit, xlimit))
	fig_title=@sprintf "ising D=%i, β=%i, freq_num=%i" D β freq_num
	plot!(xlabel="ω", ylabel="χ''(ω)/ω", title=fig_title)
end

# ╔═╡ d2e7f35f-5ef1-402c-b4ed-62682515edc9
begin
	plot(omega, de, line=(:black,2),label="theory")
	plot!(dl[t][:,1], dl[t][:,2].*dl[t][:,1], line=(:solid, 2),label="Li")
	label12 = @sprintf "χ=%i" D
	label22 = @sprintf "χ=%i×2" D
	plot!(omega, dc[t,1].* omega,line=(:dash, 2),label=label12)
	plot!(omega, dc[t,2].* omega,line=(:dash, 2),label=label22)
	plot!(xlim=(0, xlimit2), ylim=(0,maximum(de)*1.3))
	fig_title2=@sprintf "ising D=%i, β=%i, freq_num=%i" D β freq_num
	plot!(xlabel="ω", ylabel="χ''(ω)", title=fig_title2)
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
# ╟─657a6556-a153-4a62-8294-a157401cfca5
# ╟─6b130640-5f2b-49fb-8520-968f840af084
# ╟─547b3223-0cf5-4331-a5ef-0fe9b64bd16d
# ╠═8c2bd68e-1ba4-46e1-a07c-a749e1e4d1dc
# ╟─75fdb0f0-df02-4567-80f0-165c5e552d99
# ╟─a44c4474-f9ac-4b6a-a861-19c786735642
# ╟─80060631-07fd-45f6-b0e2-585e65325c0c
# ╟─02674a34-9dab-4431-bb18-6e3d7e7dc61b
# ╟─d2e7f35f-5ef1-402c-b4ed-62682515edc9
# ╟─2582d5b5-1bbd-4e68-8324-77f26af0116d
# ╟─c7972842-afb9-4440-b35f-a7ffe45fa5b7
# ╠═ded77134-54d7-4df9-b573-1830cb0b5295
# ╟─c36eb877-d015-4947-8d95-61538a636646
# ╟─7a9d8aa3-fd9d-43d0-bec6-92dcb59df920
# ╟─ec19d6ec-2a4e-468b-8b87-d45a25d9a747
# ╟─0be8ccc2-00ab-46c6-8aa9-f7173d2c89cd
# ╟─f581e951-75ca-402a-876e-83c44939737c
# ╟─8881d4a0-8bdf-4538-89a5-110d92c4f759
# ╟─3d7f5cf2-6a6e-4900-b4f4-bdfab9997dcf
# ╟─0c6f07df-9192-4d31-bb03-271fcf56864f
# ╟─7b3e2635-5a3f-4a65-b8bb-f2170155677b
# ╟─3fa261a6-5334-11ec-1dd7-8b427445a7ab
