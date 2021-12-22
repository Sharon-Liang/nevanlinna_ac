using LinearAlgebra
using Random; Random.seed!()
using HCubature
using DelimitedFiles, Printf
using Plots

function Masubara_freq(n::Int64, β::Real; type::Symbol= :f)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end


"""
G(iωn) = 1/2π ∫dΩ A(Ω)/(iωn - Ω)
"""
function Masubara_GF(n::Int64, A::Function, β::Real;
    type::Symbol = :f, Λ::Float64=100., err::Float64=eps())
    ωn = Masubara_freq(n, β, type=type)
    res = hquadrature(ω -> A(ω)/(1.0im*ωn-ω), -Λ, Λ,rtol=err)
    return res[1]/(2π)
end


"""
Gaussian approximation of Delta function
"""
function delta(x::Real; c::Real = 0., η::Real = 0.05)
    num = η
    den = ((x-c)^2 + η^2) * π
    return num/den
end


"""
    Gaussian functions: 
    f(x) = 1/(σ√(2π)) exp(-((x-μ)/σ)^2 /2 )
"""
function gaussian(x::Real, μ::Real, s2::Real)
    if s2 ≤ 0 
        @error "σ^2 should be positive!"
    else
        norm = √(2π * s2)
        m = (x-μ)^2 / s2
        return 1/norm*exp(-m/2)
    end
end

function multi_gaussian(N::Int64; μrange::Vector{Float64}=[-3.,3.], s2max::Real=3.)
    μs, s2 = rand(Float64, N), rand(Float64, N)
    μs = μs .* (μrange[2] - μrange[1]) .+ μrange[1]
    s2 = s2 .* s2max
    function sum_gaussian(x::Real)
        res = 0.
        for i = 1: N
            res += gaussian(x, μs[i], s2[i])
        end
        res/N
    end
    return sum_gaussian, μs, s2
end


"""
    Generate test data
"""
omega = [i for i in range(-π, π, length=200)]
Fn = 40
β = 20
ωn = [Masubara_freq(n,β) for n=1:Fn]

p1 = "./data/gaussian/giwn_delta.txt"
op1 = "./data/gaussian/A_delta.txt"

f = x -> delta(x, c = 1)
A = [f(ω) for ω in omega]
G = [Masubara_GF(n,f,β) for n=1:Fn]

plot(omega, A, lw=2)

open(p1, "w") do file 
    for i=1:Fn
        writedlm(file, [ωn[i] real(G[i]) imag(G[i])])
    end
end

open(op1, "w") do file 
    for i=1:length(omega)
        writedlm(file, [omega[i] A[i]])
    end
end


N = 1
p2 = @sprintf "./data/gaussian/giwn_gaussian_%i.txt" N
op2 = @sprintf "./data/gaussian/A_gaussian_%i.txt" N

f, _, _ = multi_gaussian(N)
A = [f(ω) for ω in omega]
G = [Masubara_GF(n,f,β) for n=1:Fn]

plot(omega, A, lw=2)

open(p2, "w") do file 
    for i=1:Fn
        writedlm(file, [ωn[i] real(G[i]) imag(G[i])])
    end
end

open(op2, "w") do file 
    for i=1:length(omega)
        writedlm(file, [omega[i] A[i]])
    end
end