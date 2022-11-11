using NevanlinnaAC
using NevanlinnaAC: ngradient, gradient_function
using NevanlinnaAC: hardy_expand
using NevanlinnaAC: loss_fermi, loss_bose

using Test, Zygote, Optim
using DelimitedFiles, Printf
using Random; Random.seed!()
using FiniteDifferences

fermi_clean_data = readdlm("./data/gaussian/giwn_gaussian_3.txt")
bose_clean_data = readdlm("./data/giwn_hubbard_M_point_beta_1.00.txt", skipstart = 2)

Nf = 10
fermiX = convert(Vector{ComplexF64}, fermi_clean_data[1:Nf,1])
fermiY = fermi_clean_data[1:Nf,2] + 1.0im * fermi_clean_data[1:Nf,3]

boseX = convert(Vector{ComplexF64}, bose_clean_data[2:Nf+1,1])
boseY = bose_clean_data[2:Nf+1,2] + 1.0im * bose_clean_data[2:Nf+1,3]

@testset "hardy_expand" begin
    z = rand() + 0.01im
    f(params) = imag(hardy_expand(z, params))
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

@testset "generalized_schur" begin
    x, y = toGeneralizedSchurdata(fermiX, fermiY, Fermi)
    z = rand() + 0.01im
    f(params) = generalized_schur(z, x, y, params) |> imag
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz; rtol=1.e-7)
end

@testset "nevanlinna" begin
    x, y = toNevanlinnadata(fermiX, fermiY, Fermi)
    z = rand() + 0.01im
    f(params) = nevanlinna(z, x, y, params) |> imag
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

@testset "spectral_function_value_fermi" begin
    ω = rand()
    f(params) = spectral_function_value_fermi(ω, fermiX, fermiY, params)
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

@testset "spectral_function_value_bose" begin
    ω = rand()
    g0 = 0.0
    f(params) = spectral_function_value_bose(ω, boseX, boseY, g0, params)
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

@testset "loss_fermi gradient" begin
    p0 = zeros(6)
    loss = params -> loss_fermi(params, fermiX, fermiY)
    gn = ngradient(loss, p0)[1]
    gz = gradient(loss, p0)[1]
    @test ≈(gn, gz; rtol=1.e-6)
end

@testset "loss_bose " begin
    p0 = rand(6)
    g0 = 0.0
    setprecision(BigFloat, 128)
    for isodd in [true, false], use_g0 in [true, false]
        loss = params -> loss_bose(params, boseX, boseY, g0; isodd, float_type = BigFloat, use_g0)

        gn = ngradient(loss, p0)[1]
        gz = gradient(loss, p0)[1]
        @test ≈(gn, gz, rtol=1.e-6)
    end
end