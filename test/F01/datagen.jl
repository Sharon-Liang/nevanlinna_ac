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
    delta(x::Real; c::Real = 0., η::Real = 0.05)

Gaussian approximation of Delta function
"""
function delta(x::Real; c::Real = 0., η::Real = 0.05)
    num = η
    den = ((x-c)^2 + η^2) * π
    return num/den
end


#Generate spectrum
β = 10
wmax = π
nmesh = 200
wmesh = [i for i in range(-wmax, wmax, length = nmesh)]

c = rand()
A = x->delta(x; c)

Aw = map(A, wmesh)
open("/home/sliang/JuliaCode/NevanlinnaAC/test/F01/Aw.txt", "w") do io
    writedlm(io, [wmesh Aw])
end

ngrid = 20
wn, giwn = Masubara_GF(ngrid, A, β, Fermi)
open("/home/sliang/JuliaCode/NevanlinnaAC/test/F01/giwn.txt", "w") do io
    writedlm(io, [wn real.(giwn) imag.(giwn)])
end

