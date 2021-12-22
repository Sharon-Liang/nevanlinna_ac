
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
function spectrum(ω::Real, η::Real, 
    x::Vector, y::Vector, type::Symbol;
    optim = :none)
    x, y = toNevanlinnadata(x,y,type)
    type == :f ? name = "A(ω)" : name = "ωA(ω)"
    z = ω + 1.0im * η
    res = nevanlinna(z, x, y, optim=optim)
    res = 2*res |> imag
    return res, name
end