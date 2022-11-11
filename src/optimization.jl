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
function loss_bose(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, isodd::Bool=false, kwargs...)
    if rem(Nω, 2)==0 Nω = Nω + 1 end
    ωlist, Alist = spectral_function(Bose, x, y, g0, params; ωmax, Nω, kwargs...)
    
    L = 2*ωmax; Δω = L/Nω
    Ãlist = Alist ./ ωlist
    sum_rule = abs(-g0 - sum(Ãlist) * Δω/2π)^2

    dA2_dω2 = fft_derivative(Alist, L, 2)
    smooth_condition = λ * norm(dA2_dω2)^2 * Δω

    if isodd
        odd_condition = norm((Alist[2:end] .+ reverse(Alist[2:end]))./2)^2 * Δω
        return sum_rule + smooth_condition + odd_condition
    else
        return sum_rule + smooth_condition
    end
end