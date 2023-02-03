#delta spectrum

using LinearAlgebra
using Random; Random.seed!()
using HCubature
using DelimitedFiles
using NevanlinnaAC

"""
    make_grid(ngrid::Int64, β::Real, type::OperatorType)

Make Masubara frequency grid
"""
function make_grid(ngrid::Int64, β::Real, type::OperatorType)
    if type == Fermi
        grid = [(2n+1)*π/β  for n = 0: ngrid-1]
    else
        grid = [2n*π/β  for n = 0: ngrid-1]
    end
    
    return grid  
end


"""
    Masubara_GF(n::Int64, Aω::Function, β::Real, type::OperatorType)

Calculate Masubara Green's function ``G(iωₙ) = ∫dΩ A(Ω)/(iωₙ - Ω)`` for given ``A(Ω)``.
"""
function Masubara_GF(ngrid::Int64, A::Function, β::Real, type::OperatorType; Λ::Float64=100., err::Float64=√eps())
    grid = make_grid(ngrid, β, type)
    func = ωₙ -> hquadrature(Ω -> A(Ω)/(1.0im*ωₙ-Ω), -Λ,  Λ, rtol=err)[1]
    
    return grid, map(func, grid)
end


"""
    gaussian(x::Real, μ::Real, σ²::Real)

Gaussian function: 
    f(x) = 1/(σ√(2π)) exp(-((x-μ)/σ)^2 /2 )
"""
function gaussian(x::Real, μ::Real, σ²::Real)
    if σ² ≤ 0 
        @error "σ² should be positive!"
    else
        norm = √(2π * σ²)
        m = (x-μ)^2 / σ²
        return 1.0/norm*exp(-m/2)
    end
end


function multi_gaussian(N::Int64; μrange::Vector=[-3.,3.], σ²max::Real=3.)
    μs, σ² = rand(Float64, N), rand(Float64, N)
    μs = μs .* (μrange[2] - μrange[1]) .+ μrange[1]
    σ² = σ² .* σ²max
    function sum_gaussian(x::Real)
        res = 0.
        for i = 1: N
            res += gaussian(x, μs[i], σ²[i])
        end
        res/N
    end
    return sum_gaussian, μs, σ²
end


#Generate spectrum
β = 10
wmax = 2π
nmesh = 200
wmesh = [i for i in range(-wmax, wmax, length = nmesh)]

N = rand([3,4,5])
A, _, _ = multi_gaussian(N)

Aw = map(A, wmesh)
open("/home/sliang/JuliaCode/NevanlinnaAC/test/F03/Aw.txt", "w") do io
    writedlm(io, [wmesh Aw])
end

ngrid = 20
wn, giwn = Masubara_GF(ngrid, A, β, Fermi)
open("/home/sliang/JuliaCode/NevanlinnaAC/test/F03/giwn.txt", "w") do io
    writedlm(io, [wn real.(giwn) imag.(giwn)])
end

