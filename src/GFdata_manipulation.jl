"""
    Read Masubara Green's function data
"""
function readGF(path::String; reverse::Bool=true, num::Int64 = 100)
    d = readdlm(path)
    x = 1.0im * d[:,1]
    y = d[:,2] + 1.0im * d[:,3]
    if reverse 
        x = reverse(x); y=reverse(y)
    end
    num = min(length(x), num)
    return x[1:num], y[1:num]
end

"""
    Convert Masubara frequency Green's data to Nevanlinna data
    input x: ωn
    input y: G(iωn), complex
    
    output x: iωn
    output y: -G(iωn) for fermions
              -iiωn*G(iωn) for bosons
"""
function toNevanlinnadata(x::Vector, y::Vector, type::Symbol)
    x = 1.0im * x
    if type == :f
        y = -y
    elseif type == :b
        y = -x .* y
    else
        @error "type should be :f for fermions and :b for bosons"
    end
    return x, y
end

"""
    Spectrum of Nevanlinna function G(z)
    for fermions, it is A(ω)
    for bosons, it is ωA(ω)
"""
function spectrum(ω::Vector{T} where T<:Real, η::Real, 
    x::Vector, y::Vector, type::Symbol;
    optim = :none)
    x, y = toNevanlinnadata(x,y,type)
    if isNevanlinnasolvable(x,y)[1] == false @warn "Nevanlinna unsolvable!" end
    type == :f ? name = "A(ω)" : name = "ωA(ω)"
    z = ω .+ 1.0im * η
    res = [nevanlinna(i, x, y, optim=optim) for i in z]
    res = imag.(2*res)
    return res, name
end

function spectrum(ω::Real, η::Real, 
    x::Vector, y::Vector, type::Symbol;
    optim = :none)
    return spectrum([ω], η, x, y, type, optim = optim)
end

