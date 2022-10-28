Zygote.@nograd ispossemidef
Zygote.@nograd pick_matrix
Zygote.@nograd isGeneralizedSchursovable
Zygote.@nograd isNevanlinnasolvable
Zygote.@nograd schur_parameter
Zygote.@nograd coefficient

Zygote.@nograd fftfreq

#TODO: fermionic case 
function loss_fermi(params::AbstractArray, x::AbstractArray, y::AbstractVector; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)
    _, Alist = spectral_function(Fermi, x, y, params; ωmax, Nω, kwargs...)

    L = 2*ωmax
    klist = fftfreq(Nω)*Nω
    fft_Alist = fft(convert(Vector{Float64}, Alist))
    dA2_dω2 = real.(ifft((2π*im/L .* klist).^2 .* fft_Alist));  

    loss = abs(1.0 - sum(Alist)*L/Nω)^2 + λ * norm(dA2_dω2)^2
    return loss
end

