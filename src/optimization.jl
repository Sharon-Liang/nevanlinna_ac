Zygote.@nograd pick_matrix
Zygote.@nograd issolvable
Zygote.@nograd schur_parameter
Zygote.@nograd coefficient

Zygote.@nograd toNevData
Zygote.@nograd toGenSchurData
Zygote.@nograd fftfreq

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
    fft_derivative(d::AbstractVector{<:Real}, L::Real, order::Real=1)

    calculate the `n`-th order derivative of given data ``d``. ``L`` is the lenth of definition range `x ∈ [xmin, xmin+L)`
"""
#TODO:Fail, read FFT and rewrite this part
#TODO: modify make_mesh function accordingly
function fft_derivative(d::AbstractVector{<:Real}, L::Real, order::Real=1)
    nmesh = lastindex(d)
    klist = fftfreq(nmesh) * nmesh
    if eltype(datalist) <: FFTW.fftwNumber
        fft_datalist = fft(datalist)
    else
        fft_datalist = fft(convert(Vector{Float64}, datalist))
    end
    return real.(ifft((2π*im/L .* klist).^n .* fft_datalist))
end


function loss(params::AbstractArray, d::RawData, option::Options; λ::Real=1.e-4)
    @unpack wmax, otype = option

    wmesh, Aw = spectral_function(option, d, params)
    L = 2*wmax
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

