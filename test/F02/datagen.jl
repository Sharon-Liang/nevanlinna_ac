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
function Masubara_GF(ngrid::Int64, A::Function, β::Real, type::OperatorType; Λ::Float64=100., err::Float64=1.e-15)
    grid = make_grid(ngrid, β, type)
    func = ωₙ -> hquadrature(Ω -> A(Ω)/(1.0im*ωₙ-Ω), -Λ,  Λ, rtol=err)[1]
    
    return grid, map(func, grid)
end


"""
    gaussian(x::Real, μ::Real, s2::Real)

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


#Generate spectrum
β = 10
wmax = π
nmesh = 200
wmesh = [i for i in range(-wmax, wmax, length = nmesh)]

μ = rand()
σ² = rand()

A = x->gaussian(x, μ, σ²)

Aw = map(A, wmesh)
open("/home/sliang/JuliaCode/NevanlinnaAC/test/F02/Aw.txt", "w") do io
    writedlm(io, [wmesh Aw])
end

ngrid = 20
wn, giwn = Masubara_GF(ngrid, A, β, Fermi)
open("/home/sliang/JuliaCode/NevanlinnaAC/test/F02/giwn.txt", "w") do io
    writedlm(io, [wn real.(giwn) imag.(giwn)])
end

