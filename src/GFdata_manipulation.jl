"""
    Read Masubara Green's function data
"""
function readGF(ftype::DataType, path::String; rev::Bool=true, num::Integer = 100)
    Im = one(ftype)im  #define imaginary unit
    ctype =  Complex{ftype} 
    d = readdlm(path)
    x = ctype.(d[:,1]) 
    y = ftype.(d[:,2]) + Im * ftype.(d[:,3])
    num = min(length(x), num)
    x1, y1 = x[1:num], y[1:num]
    if rev == true
        x1 = reverse(x1); y1=reverse(y1)
    end
    return x1, y1
end

function readGF(path::String; rev::Bool=true, num::Integer = 100)
    readGF(Float64, path; rev = rev, num = num)
end


"""
    Convert Masubara frequency Green's data to Nevanlinna data
    input x: ωn
    input y: G(iωn), complex
    
    output x: iωn
    output y: -G(iωn) for fermions
              -iiωn*G(iωn) for bosons
"""
function toNevanlinnadata(ftype::DataType, x::AbstractVector{T}, y::AbstractVector{T}, operator_type::OperatorType) where T
    ctype = Complex{ftype}
    Im = one(ftype)im
    x = ctype.(x); y = ctype.(y)
    x = Im * x
    if operator_type == Fermi
        y = -y
    else
        y = -x .* y
    end
    return x, y
end

function toNevanlinnadata(x::AbstractVector{T}, y::AbstractVector{T}, operator_type::OperatorType) where T
    Im = one(T)im
    x = Im * x
    if operator_type == Fermi
        y = -y
    else
        y = -x .* y
    end
    return x, y
end


"""
    Convert Masubara frequency Green's data to Generalized Schur data
    input x: ωn
    input y: G(iωn), complex
    
    output x: iωn
    output y: mt(-G(iωn), 1im) for fermions
              mt(-iiωn*G(iωn, 1im)) for bosons
"""
function toGeneralizedSchurdata(ftype::DataType, x::AbstractVector{T}, y::AbstractVector{T}, operator_type::OperatorType) where T
    x, y = toNevanlinnadata(ftype, x, y, type)
    return x, mti.(ftype, y)
end

function toGeneralizedSchurdata(x::AbstractVector{T}, y::AbstractVector{T}, operator_type::OperatorType) where T
    x, y = toNevanlinnadata(x, y, type)
    return x, mti.(y)
end


"""
    Spectrum of Nevanlinna function G(z)
    for fermions, it is A(ω)
    for bosons, it is ωA(ω)
"""
function spectrum(ω::Vector{T} where T<:Real, η::Real, 
    x::AbstractVector, y::AbstractVector, operator_type::OperatorType;
    optim = :none)
    x, y = toNevanlinnadata(x,y,operator_type )
    if isNevanlinnasolvable(x,y)[1] == false @warn "Nevanlinna unsolvable!" end
    operator_type == Fermi ? name = "A(ω)" : name = "ωA(ω)"
    z = ω .+ 1.0im * η
    res = zeros(eltype(y), length(ω))
    for i = 1:length(ω)
        res[i] = nevanlinna(z[i], x, y, optim=optim) 
    end
    res = imag.(2*res)
    return res, name
end

function spectrum(ω::Real, η::Real, 
    x::AbstractVector, y::AbstractVector, operator_type::OperatorType;
    optim = :none)
    return spectrum([ω], η, x, y, operator_type, optim = optim)
end

