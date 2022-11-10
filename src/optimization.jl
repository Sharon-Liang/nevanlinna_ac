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
    
    L = 2*ωmax; Δω = L/Nω
    dA2_dω2 = fft_derivative(Alist, L, 2)

    sum_rule = abs(1.0 - sum(Alist) * Δω)^2
    smooth_condition = λ * norm(dA2_dω2)^2 * Δω
 
    return sum_rule + smooth_condition
end


"""
loss_bose(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)

"""
function loss_bose_fail(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, isodd::Bool=true, kwargs...)
    A_function = ω -> spectral_function_value_bose(ω, x, y, g0, params; kwargs...)/2π
    rho_function = ω -> A_function(ω)/ω
    dA_function = ω -> ngradient(A_function, ω)

    int_rho_positive =  hquadrature(rho_function, 0, ωmax)[1]
    int_rho_negative = hquadrature(rho_function, -ωmax, 0)[1]
    sum_rule = abs(-g0 - (int_rho_positive + int_rho_negative))^2

    smooth_condition = λ * hquadrature(dA_function, -ωmax, ωmax)[1]

    if isodd
        odd_condition = hquadrature(ω-> norm(0.5 * (A_function(ω) + A_function(-ω))), -ωmax, ωmax)[1]
        return  sum_rule + smooth_condition + odd_condition
    else
        return sum_rule + smooth_condition
    end
end


function loss_bose(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, isodd::Bool=true, kwargs...)
    ωlist, Alist = spectral_function(Bose, x, y, g0, params; ωmax, Nω, kwargs...)

    L = 2*ωmax; Δω = L/Nω
    klist = fftfreq(Nω)*Nω
    fft_Alist = fft(convert(Vector{Float64}, Alist))
    dA_dω = real.(ifft((2π*im/L .* klist) .* fft_Alist))
    dA2_dω2 = real.(ifft((2π*im/L .* klist).^2 .* fft_Alist))

    smooth_condition = λ * norm(dA2_dω2)^2 * Δω

    ind = findall(x -> x == 0., ωlist)[1]
    s0 = dA_dω[ind]
    Ãlist = vcat(Alist[1:ind-1] ./ ωlist[1:ind-1], [s0], Alist[ind+1:end] ./ ωlist[ind+1:end])
    sum_rule = abs(-g0 - sum(Ãlist) * L/Nω /2π)^2

    return sum_rule + smooth_condition
end