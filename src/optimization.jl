Zygote.@nograd pick_matrix
Zygote.@nograd issolvable
Zygote.@nograd schur_parameter
Zygote.@nograd coefficient

Zygote.@nograd toNevData
Zygote.@nograd toGenSchurData
Zygote.@nograd fftfreq
Zygote.@nograd make_mesh

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
    N = lastindex(d)
    k = fftfreq(N) * N * (2π*im/L)
    if eltype(d) <: FFTW.fftwNumber
        fft_d = fft(d)
    else
        fft_d = fft(convert(Vector{Float64}, d))
    end

    res = ifft(k.^order .* fft_d)
    return real.(res)
end


function loss(params::AbstractArray, d::RawData, option::Options; λ::Real=1.e-4)
    @unpack wmax, otype = option

    wmesh, Aw = spectral_function(option, d, params)
    L = 2*wmax
    Δω = L / lastindex(wmesh)

    ∂²Aw = fft_derivative(Aw, L, 2)
    smooth_condition = λ * norm(∂²Aw)^2 * Δω 

    #smooth_condition = λ * norm(Aw)
    
    if otype == Fermi
        sum_rule = abs(1.0 - sum(Aw) * Δω)^2
    else
        ind = findall(x -> x==0.0, grid)[1]
        Ãw = @. Aw / wmesh
        sum_rule = abs(-d.value[ind] - sum(Ãw) * Δω)^2
    end

    return smooth_condition + sum_rule
end

