using nevanlinna_ac
import nevanlinna_ac._mti
include("optim")
using LinearAlgebra, Optim, Zygote, FFTW
using ChainRulesCore
using Plots

"""
    H2basis(z, k)

The `k`-th orthonormal basis of `H2` space.
"""
function H2basis(z::Number, k::Int64)
    return 1/(√π*(z+oneunit(z)im)) * _mti(z)^k
end


"""
    hardy_expand(z, Nh::Int64,[, ak::AbstractVector = rand(typeof(z), Nh), bk::AbstractVector = rand(typeof(z), Nh)])

    generate a function `f(z)` by the Hardy basis up to order ``Nk``, the corresponding coefficients are ``ak`` and ``bk``.
"""
function hardy_expand(z::Number, ak::AbstractVector, bk::AbstractVector)
    if length(ak) != length(bk) @error "Length of coefficient Vectors should be the same!" end
    Hklist = [H2basis(z, k) for k=1:length(ak)]
    res = sum(ak .* Hklist .+ bk .* conj.(Hklist))
    return res
end


"""
    Loss function 
"""
function _loss_fermi(params::AbstractArray;  x::AbstractVector, y::AbstractVector, float_type::DataType=Double64, η::Real = 0.05, toreverse::Bool=true, init_func::Function)
    init_func = z -> init_func(z, params)
    ωlist, Alist = spectral_function(Fermi, x, y; init_func)
    
    L1 = abs(1.0 - sum(Alist) * L/Nω)^2

    klist = fftfreq(Nω)*Nω
    fft_Alist = fft(convert(Vector{Float64}, Alist))
    dA2_dw2 = real.(ifft((2π*im/L .* klist).^2 .* fft_Alist));

    L2 = λ * norm(dA2_dw2)

    return L1 + L2
end











"""
    spectral_function_hardy_optim(operator_type::OperatorType, x::AbstractVector, y::AbstractVector[, float_type::DataType=Double64]; ωmax::Real = 4π, Nω::Integer = 500, η::Real = 1.e-5, toreverse::Bool=true, hardy_order::Int64 = 25)

Calculate the spectral function `A(ω)` for given dataset `{x=ωn, y=G(iωn)}` in the range `[-ωmax, ωmax)` with `Nω` discrete values.
"""
#function spectral_function_hardy_optim(operator_type::OperatorType, x::AbstractVector, y::AbstractVector, float_type::DataType=Double64; ωmax::Real = 4π, Nω::Integer = 500, η::Real = 1.e-5, toreverse::Bool=true, hardy_order::Int64 = 25, λ::Real = 1.e-4)
    ωmax = 4π; Nω = 500
    L = 2*ωmax
    float_type = Float64
    hardy_order = 25
    operator_type = Fermi
    η = 1.e-5
    toreverse = false

    using DelimitedFiles
    A1data = readdlm("/Users/liangshuang/Desktop/CAS-code/nevanlinna_ac/data/gaussian/A_gaussian_1.txt")
    ωlist = A1data[:,1]
    A1 = A1data[:, 2]
    xdata, ydata = readGF("/Users/liangshuang/Desktop/CAS-code/nevanlinna_ac/data/gaussian/giwn_gaussian_1.txt");
    λ = 1.e-4

    wlist, A1_ac = spectral_function(Fermi, xdata, ydata, float_type, toreverse = false);

    p1 = plot(ωlist, A1, line=2, label="exact")
    plot!(wlist, A1_ac, line=(1,:dash), marker=2, label="ac")

    ωlist = convert(Vector{float_type}, (-Nω/2:Nω/2-1)*L/Nω)
    #define a hardy basis hardy basis 
    ak = rand(float_type, hardy_order)
    bk = rand(float_type, hardy_order)

    loss1 = x0 -> begin
        ak = x0[:,1]
        bk = x0[:,2]
        init_func = z -> hardy_expand(z, hardy_order, ak, bk)
        _, Alist = spectral_function(Fermi, xdata, ydata, float_type; init_func,toreverse = false);
        klist = fftfreq(Nω)*Nω
        fft_Alist = fft(convert(Vector{Float64}, Alist))
        dA2_dw2 = real.(ifft((2π*im/L .* klist).^2 .* fft_Alist));  
        res =  abs(1.0 - sum(Alist) * L/Nω)^2 + 1.e-4 * norm(dA2_dw2)
        return res
    end

    result = optimize(loss1, hcat(ak,bk), show_trace = true)

    init_func = z -> hardy_expand(z, hardy_order, result.minimizer[:,1], result.minimizer[:,2])
    A1_ac_hardy = map(ω -> spectral_function_value(ω, operator_type, xdata, ydata, float_type; η, toreverse, init_func), ωlist)
    plot!(wlist, A1_ac_hardy, line=(1,:dash), marker=2, label="ac_hardy")











#module OptimFunctions
#Reference: https://github.com/baggepinnen/FluxOptTools.jl/blob/master/src/FluxOptTools.jl
veclength(grads::Zygote.Grads) = sum(length(grads[p]) for p in grads.params)
veclength(params::Zygote.Params) = sum(length, params.params)
veclength(x) = length(x)

Base.zeros(grads::Zygote.Grads) = zeros(veclength(grads))
Base.zeros(pars::Zygote.Params) = zeros(veclength(pars))


"""
    optim_functions(loss, pars::Zygote.Params)

Generate the loss function, gradient function and and `p0`, a vectorized version of pars.
"""
function optim_functions(loss, pars::Zygote.Params)
    grads = Zygote.gradient(loss, pars)
    p0 = zeros(pars)
    copy!(p0, pars)

    gradient_function = function (g,w)
        copy!(pars, w)
        grads = Zygote.gradient(loss, pars)
        copy!(g, grads)
    end

    loss_function = function (w)
        copy!(pars, w)
        loss()
    end
    return p0, loss_function, gradient_function
end








