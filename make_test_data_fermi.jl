using LinearAlgebra
using Random; Random.seed!()
using HCubature
using DelimitedFiles, Printf
using Plots
using NevanlinnaAC: OperatorType, Bose, Fermi

function Masubara_freq(n::Int64, β::Real, type::OperatorType)
    type == Bose ? N = 2n : N = 2n + 1
    return N*π/β
end


"""
G(iωn) = 1/2π ∫dΩ A(Ω)/(iωn - Ω)
"""
function Masubara_GF(n::Int64, A::Function, β::Real, type::OperatorType; Λ::Float64=100., err::Float64=1.e-15)
    ωn = Masubara_freq(n, β, type)
    res = hquadrature(ω -> A(ω)/(1.0im*ωn-ω), -Λ,  Λ, rtol=err)
    return res[1]/2π
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

function multi_gaussian(N::Int64; μrange::Vector=[-3.,3.], s2max::Real=3.)
    μs, s2 = rand(Float64, N), rand(Float64, N)
    μs = μs .* (μrange[2] - μrange[1]) .+ μrange[1]
    s2 = s2 .* s2max
    function sum_gaussian(x::Real)
        res = 0.
        for i = 1: N
            res += gaussian(x, μs[i], s2[i]) * 2π
        end
        res/N
    end
    return sum_gaussian, μs, s2
end


"""
    Generate test data
"""
omega = [i for i in range(-4π, 4π, length=500)]
Fn = 40
β = 20
ωn = [Masubara_freq(n, β, Bose) for n=1:Fn]




f, _, _ = multi_gaussian(3, μrange=[1,3])






p1 = "./data/gaussian/giwn_delta_eta_0.05.txt"
op1 = "./data/gaussian/A_delta_eta_0.05.txt"

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


N = 3
p2 = @sprintf "./data/gaussian/giwn_gaussian_%i.txt" N
op2 = @sprintf "./data/gaussian/A_gaussian_%i.txt" N

f, _, _ = multi_gaussian(N, μrange=[-3,3])
A = [f(ω) for ω in omega]
plot(omega, A, lw=2)

G = [Masubara_GF(n,f,β, err=1.e-15) for n=1:Fn]



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

println("finish!")