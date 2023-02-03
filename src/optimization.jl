Zygote.@nograd pick_matrix
Zygote.@nograd issolvable
Zygote.@nograd schur_parameter
Zygote.@nograd coefficient

Zygote.@nograd toNevData
Zygote.@nograd toGenSchurData
Zygote.@nograd fftfreq

"""
    ngradient(f, xs::AbstractArray...)

Calculate the gradient of `f(xs...)` by finite differences.
"""
function ngradient(f, xs::AbstractArray...; step::Real=sqrt(eps()))
    #https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
    grads = zero.(xs)
    for (x, Δ) in zip(xs, grads), i in eachindex(x)
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
        fft_datalist = fft(datalist)
    else
        fft_datalist = fft(convert(Vector{Float64}, datalist))
    end
    return real.(ifft((2π*im/L .* klist).^n .* fft_datalist))
end


function loss(params::AbstractArray, d::RawData, option::Options; λ::Real=1.e-4)
    @unpack wmax, wmin, otype = option

    wmesh, Aw = spectral_function(option, d, params)
    L = wmax - wmin
    Δω = L / lastindex(wmesh)

    ∂²Aw = fft_derivative(Aw, L, 2)
    smooth_condition = λ * norm(∂²Aw)^2 * Δω 
    
    if otype == Fermi
        sum_rule = abs(1.0 - sum(Aw) * Δω)^2
    else
        ind = findall(x -> x==0.0, grid)
        Ãw = @. Aw / wmesh
        sum_rule = abs(-d.value[ind] - sum(Ãw) * Δω)^2
    end

    return smooth_condition + sum_rule
end

