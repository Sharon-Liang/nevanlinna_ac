using DelimitedFiles, Printf
using Test, Random; Random.seed!()

using NevanlinnaAC
using NevanlinnaAC: ngradient, ngradient2, gradient_function, fft_derivative

using LinearAlgebra, FiniteDifferences, HCubature

bose_clean_data = readdlm("./data/giwn_hubbard_M_point.txt", skipstart = 2)

Nf = 10
boseX = convert(Vector{ComplexF64}, bose_clean_data[2:Nf+1,1])
boseY = bose_clean_data[2:Nf+1,2] + 1.0im * bose_clean_data[2:Nf+1,3]
g0 = bose_clean_data[1,2]

prec = 128
setprecision(BigFloat, prec); float_type = BigFloat


@testset "FFT and FiniteDifferences - 1st order gradient" begin
    Nω = 1000
    L = 4π
    ωlist = convert(Vector{float_type}, (-Nω/2:Nω/2-1)*L/Nω)

    A_function = ω -> spectral_function_value_bose(ω, boseX, boseY, g0; use_g0 =false, float_type)

    gn1 = map(ω -> ngradient(A_function, ω), ωlist)
    gn2 = fft_derivative(map(A_function, ωlist), L, 1)
    gn3 = map(ω -> central_fdm(5,1)(A_function, ω), ωlist)

    @test gn1 ≈ gn2 rtol=0.01
    @test gn1 ≈ gn3 rtol=1.e-8
end

@testset "FFT and FiniteDifferences - 2nd order gradient" begin
    Nω = 1000
    L = 4π
    ωlist = convert(Vector{float_type}, (-Nω/2:Nω/2-1)*L/Nω)

    A_function = ω -> spectral_function_value_bose(ω, boseX, boseY, g0; use_g0 =false, float_type)

    gn1 = map(ω -> ngradient2(A_function, ω), ωlist)
    gn2 = fft_derivative(map(A_function, ωlist), L, 2)
    gn3 = map(ω -> central_fdm(7,2)(A_function, ω), ωlist)

    @test gn1 ≈ gn2 rtol=0.02
    @test gn1 ≈ gn3 rtol=0.02 
    @test gn2 ≈ gn3 rtol=0.02
end


@testset "sum rule - HCubature" begin
    A_function = ω -> spectral_function_value_bose(ω, boseX, boseY, g0, rand(10); use_g0 = false, float_type)

    Nfunction = ω -> ngradient(A_function, ω)
    A_sum, err = hquadrature(Nfunction, -4π, 4π, atol=1.e-5)
    A_sum1, err = hquadrature(Nfunction, -4π, 0, atol=1.e-5)
    A_sum2, err = hquadrature(Nfunction, 0, 4π, atol=1.e-5)
    @show A_sum
    @show A_sum1 + A_sum2
    @show g0
end





