Zygote.@nograd ispossemidef
Zygote.@nograd pick_matrix
Zygote.@nograd isGeneralizedSchursovable
Zygote.@nograd isNevanlinnasolvable
Zygote.@nograd schur_parameter
Zygote.@nograd coefficient

Zygote.@nograd toNevanlinnadata
Zygote.@nograd toGeneralizedSchurdata
Zygote.@nograd fftfreq


"""
loss_fermi(params::AbstractArray, x::AbstractArray, y::AbstractVector; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)

    loss function for fermionic case
"""
function loss_fermi(params::AbstractArray, x::AbstractArray, y::AbstractVector; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)
    _, Alist = spectral_function(Fermi, x, y, params; ωmax, Nω, kwargs...)
    L = 2*ωmax

    klist = fftfreq(Nω)*Nω
    fft_Alist = fft(convert(Vector{Float64}, Alist))
    dA2_dω2 = real.(ifft((2π*im/L .* klist).^2 .* fft_Alist));  

    sum_rule = abs(1.0 - sum(Alist)*L/Nω)^2
    smooth_condition = λ * norm(dA2_dω2)^2
 
    return sum_rule + smooth_condition
end


"""
loss_bose(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)

"""
function loss_bose(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)
    A_function = ω -> spectral_function_value_bose(ω, x, y, g0, params; ωmax, Nω, kwargs...)

    ωlist = [i for i in range(-ωmax, ωmax, length = Nω)]
    Alist = map(A_function, ωlist)

    Ãlist = Alist ./ ωlist
    dA2_dω2 = map(ω->central_fdm(5,2)(A_function, ω), ωlist)

    smooth_condition = λ * norm(dA2_dω2)^2

    sum_rule = abs(-g0 - sum(Ãlist) * L/Nω /2π)^2

    odd_condition = (Alist .+ reverse(Alist))/2 |> norm

    return sum_rule + smooth_condition + odd_condition
end


function loss_bose_old(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)
    ωlist, Alist = spectral_function(Bose, x, y, g0, params; ωmax, Nω, kwargs...)

    L = 2*ωmax
    klist = fftfreq(Nω)*Nω
    fft_Alist = fft(convert(Vector{Float64}, Alist))
    dA_dω = real.(ifft((2π*im/L .* klist) .* fft_Alist))
    dA2_dω2 = real.(ifft((2π*im/L .* klist).^2 .* fft_Alist))

    smooth_condition = λ * norm(dA2_dω2)^2

    ind = findall(x -> x == 0., ωlist)[1]
    s0 = dA_dω[ind]
    Ãlist = vcat(Alist[1:ind-1] ./ ωlist[1:ind-1], [s0], Alist[ind+1:end] ./ ωlist[ind+1:end])
    sum_rule = abs(-g0 - sum(Ãlist) * L/Nω /2π)^2

    return sum_rule + smooth_condition
end