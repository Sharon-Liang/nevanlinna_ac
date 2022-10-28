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

    loss1 = abs(1.0 - sum(Alist)*L/Nω)^2
    loss2 = norm(dA2_dω2)^2
    loss =  loss1 + λ * loss2
    return loss
end

#TODO: loss bose: 奇函数限制
function loss_bose(params::AbstractArray, x::AbstractArray, y::AbstractVector, g0::Real; ωmax::Real=4π, Nω::Int=500, λ::Real=1.e-4, kwargs...)
    ωlist, Alist = spectral_function(Bose, x, y, g0, params; ωmax, Nω, kwargs...)
    L = 2*ωmax

    #sum rule
    ind = Int(Nω/2 + 1)
    sum_A  = (Alist[ind+1] - Alist[ind-1])/2
    sum_A += sum(Alist[1:ind-1] ./ ωlist[1:ind-1]) * L/Nω
    sum_A += sum(Alist[ind+1: end] ./ ωlist[ind+1: end]) * L/Nω
    loss1 = abs(-g0 - sum_A/2π)^2

    #smooth
    klist = fftfreq(Nω)*Nω
    fft_Alist = fft(convert(Vector{Float64}, Alist))
    dA2_dω2 = real.(ifft((2π*im/L .* klist).^2 .* fft_Alist));  
    loss2 = norm(dA2_dω2)^2

    #odd function
    loss3 = sum(abs.(Alist[2:ind-1] .+ Alist[ind+1: end])) + abs(Alist[ind])

    loss = loss1 + λ * loss2 + loss3
    return loss
end