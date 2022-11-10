"""
    OperatorType Bose Fermi

Operator type of physical observables
"""
@enum OperatorType Bose Fermi

"""
    ngradient(f, xs::AbstractArray...)

Calculate the gradient of `f(xs...)` by finite differences.
"""
function ngradient(f, xs::AbstractArray...; step::Real=sqrt(eps()))
    #https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
    grads = zero.(xs)
    for (x, Δ) in zip(xs, grads), i in 1:length(x)
      δ = step
      tmp = x[i]
      x[i] = tmp - δ/2
      y1 = f(xs...)
      x[i] = tmp + δ/2
      y2 = f(xs...)
      x[i] = tmp
      Δ[i] = (y2-y1)/δ
    end
    return grads
end

function ngradient(f, xs::Number; step::Real=sqrt(eps()))
    #https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
    δ = step
    y1 = f(xs - δ/2)
    y2 = f(xs + δ/2)
    return (y2-y1)/δ
end


"""
    ngradient2(f, xs::AbstractArray...)

Calculate the second order derivative of `f(xs...)` by finite differences.
"""
function ngradient2(f, xs::Number; step::Real=sqrt(eps()))
    δ = step
    df = x -> ngradient(f, x; step)
    y1 = df(xs - δ/2)
    y2 = df(xs + δ/2)
    return (y2-y1)/δ
end


"""
    gradient_function(loss, pars::AbstractArray)

Generate the gradient function of `loss`
"""
function gradient_function(loss, pars::AbstractArray)
    g! = function (g, pars)
        grads = Zygote.gradient(loss, pars)[1]
        copy!(g, grads)
    end
    return g!
end


"""
    fft_derivative(datalist::AbstractVector{<:Real}, L::Real, n::Real=1)

    calculate the `n`-th order derivative of given ``datalist``. ``L`` is the lenth of definition range `x ∈ [xmin, xmin+L]`
"""

function fft_derivative(datalist::AbstractVector{<:Real}, L::Real, n::Real=1)
    Nω = length(datalist)
    klist = fftfreq(Nω) * Nω
    if eltype(datalist) <: FFTW.fftwNumber
        fft_datalist = fft(dataliat)
    else
        fft_datalist = fft(convert(Vector{Float64}, datalist))
    end
    return real.(ifft((2π*im/L .* klist).^n .* fft_datalist))
end