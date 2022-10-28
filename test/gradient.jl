using NevanlinnaAC
using NevanlinnaAC: ngradient
using NevanlinnaAC: hardy_expand
using NevanlinnaAC: loss_fermi

using Test, Zygote, Optim
using DelimitedFiles

@testset "hardy_expand" begin
    f(params) = imag(hardy_expand(1.0 + 0.01im, params))
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

data  = readdlm("./data/gaussian/giwn_gaussian_3.txt")
N = 10
xdata = convert(Vector{ComplexF64}, data[1:N,1])
ydata = data[1:N,2] + 1.0im * data[1:N,3];

@testset "generalized_schur" begin
    x, y = toGeneralizedSchurdata(xdata, ydata, Fermi)
    z = rand() + 0.01im
    f(params) = generalized_schur(z, x, y, params) |> imag
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz; rtol=1.e-7)
end

@testset "nevanlinna" begin
    x, y = toNevanlinnadata(xdata, ydata, Fermi)
    z = rand() + 0.01im
    f(params) = nevanlinna(z, x, y, params) |> imag
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

@testset "spectral_function_value_fermi" begin
    ω = rand()
    f(params) = spectral_function_value_fermi(ω, xdata, ydata, params)
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

@testset "spectral_function_value_bose" begin
    ω = rand()
    g0 = -2
    f(params) = spectral_function_value_bose(ω, xdata, ydata, g0, params)
    p0 = rand(6)
    gn = ngradient(f, p0)[1]
    gz = gradient(f, p0)[1]
    @test ≈(gn, gz)
end

@testset "loss_fermi gradient" begin
    p0 = rand(6)
    loss = params -> loss_fermi(params, xdata, ydata)
    gn = ngradient(loss, p0)[1]
    gz = gradient(loss, p0)[1]
    @test ≈(gn, gz; rtol=1.e-6)
end

@testset "loss_fermi " begin
    @show p0 = rand(6)
    loss = params -> loss_fermi(params, xdata, ydata; Nω = 2000)
    g! = function (g, p0)
        grads = Zygote.gradient(loss, p0)[1]
        copy!(g, grads)
    end
    @show res = Optim.optimize(loss, g!, p0, LBFGS(), Optim.Options(show_trace=true))
end