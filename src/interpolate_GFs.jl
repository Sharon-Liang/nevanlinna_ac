"""
    toNevanlinnadata(x::AbstractVector, y::AbstractVector, type::OperatorType[, float_type::DataType = Double64] )

Convert Masubara frequency correlation function data to Nevanlinna data
    input x: ωn
    input y: G(iωn), complex numbers
    
    output x: z = iωn
    output y: -G(z) = -G(iωn) for fermions
              -zG(z) = (-)iωn * G(iωn) for bosons
"""
function toNevanlinnadata(x::AbstractVector, y::AbstractVector, operator_type::OperatorType, float_type::DataType = Double64)
    complex_type = Complex{float_type}
    Im = one(float_type)im
    x = convert(Vector{complex_type}, x)
    y = convert(Vector{complex_type}, y)
    x = Im * x
    if operator_type == Fermi
        y = -y
    else
        y = -x .* y
    end
    return x, y
end



"""
    toGeneralizedSchurdata(x::AbstractVector, y::AbstractVector, type::OperatorType[, float_type::DataType = Double64] )
    
Convert Masubara frequency Green's data to Generalized Schur data
    input x: ωn
    input y: G(iωn), complex numbers

    output x: z = iωn
    output y: mti[-G(z)] = mti[-G(iωn)] for fermions
              mti[-zG(z)] = mti[(-)iωn * G(iωn)] for bosons
"""
function toGeneralizedSchurdata(args...; kwargs...)
    x, y = toNevanlinnadata(args...; kwargs...)
    return x, map(_mti, y)
end


"""
    spectral_function_value_fermi(ω::Number, x::AbstractVector, y::AbstractVector[; float_type::DataType=Double64, η::Real = 0.05, toreverse::Bool=true, init_func::Function= z -> 0.0])

Calculate the spectral function `A(ω)` for given dataset `{x=ωn, y=G(iωn)}` at `ω` if `G(iωn)` is fermionic.
"""
function spectral_function_value_fermi(ω::Number, x::AbstractVector, y::AbstractVector, args...; float_type::DataType=Double64, η::Real = 0.05, toreverse::Bool=true, show_warning::Bool = false)
    x, y = toNevanlinnadata(x, y, Fermi, float_type)
    if toreverse == true
        x = reverse(x)
        y = reverse(y)
    end

    #if isNevanlinnasolvable(x,y)[1] == false @warn "Nevanlinna unsolvable!" end
    
    ω = convert(float_type, ω)
    z = ω + one(float_type)im * η
    Gω =  nevanlinna(z, x, y, args...; show_warning)

    return 2*imag(Gω)
end


"""
    spectral_function_value_bose(ω::Number, x::AbstractVector, y::AbstractVector, g0::Real[; float_type::DataType=Double64, η::Real = 0.05, toreverse::Bool=true, init_func::Function= z -> 0.0])

Calculate the spectral function `A(ω)` for given dataset `{x=ωn, y=G(iωn)}` at `ω` if `G(iωn)` is bosonic.
"""
function spectral_function_value_bose(ω::Number, x::AbstractVector, y::AbstractVector, g0::Real, args...; float_type::DataType=Double64, η::Real = 0.05, toreverse::Bool=true, alpha::Real=1.0, use_g0::Bool=false, show_warning::Bool=false)
    x, y = toNevanlinnadata(x, y, Bose, float_type)
    if toreverse == true
        x = reverse(x)
        y = reverse(y)
    end
    
    if use_g0
        g0 = convert(float_type, g0)
        x_new = vcat(x, [one(float_type)im * η])
        y_new = vcat(y, [-one(float_type)im * η * (one(float_type) - alpha * η^2) * g0])
    else
        x_new = x
        y_new = y
    end
    #if isNevanlinnasolvable(x,y)[1] == false @warn "Nevanlinna unsolvable!" end
    
    ω = convert(float_type, ω)
    z = ω + one(float_type)im * η
    Gω =  nevanlinna(z, x_new, y_new, args...; show_warning)
    
    num = ω*imag(Gω) - η*real(Gω)
    den = ω^2 + η^2
    return 2*num/den
end


"""
    spectral_function_value(operator_type::OperatorType, args...; kwargs...)

Calculate the spectral function `A(ω)` for given dataset `{x=ωn, y=G(iωn)}` at `ω`.
"""
function spectral_function_value(operator_type::OperatorType, args...; kwargs...)
    if operator_type == Bose
        return spectral_function_value_bose(args...; kwargs...)
    else
        return spectral_function_value_fermi(args...; kwargs...)
    end
end


"""
    spectral_function(operator_type::OperatorType, args...[; float_type::DataType=Double64, ωmax::Real = 4π, Nω::Integer = 500, kwargs...])

Calculate the spectral function `A(ω)` for given dataset `{x=ωn, y=G(iωn)}` in the range `[-ωmax, ωmax)` with `Nω` discrete values.
"""
function spectral_function(operator_type::OperatorType, args...; float_type::DataType=Double64, ωmax::Real = 4π, Nω::Integer = 500, kwargs...)
    L = 2*ωmax
    ωlist = convert(Vector{float_type}, (-Nω/2:Nω/2-1)*L/Nω)
    Alist = map(ω -> spectral_function_value(operator_type, ω, args...; float_type, kwargs...), ωlist)
    return ωlist, Alist
end