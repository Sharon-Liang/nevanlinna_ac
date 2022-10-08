"""
    Identity matrix
"""
function eye(n::Int64)
    return Matrix{Float64}(I,n,n)
end

function eye(dtype::DataType, n::Int64)
    return Matrix{dtype}(I,n,n)
end


"""
    linear_fractional_transform(z, Y)
    lft(z, Y)
    
Linear fractional transform transforms the close unit disk `D` in the complex plane to itself, in which point `Y` is mapped to the center of `D`. Moreover, it is a one to one mapping of the boundary of `D`.
"""
function linear_fractional_transform(z, Y) 
    #if abs(z) ≥ 1 @warn "|z|<1 is required." end
    #if abs(Y) ≥ 1 @warn "New center |Y|<1 is required." end
    dtype = eltype(z)
    num = z - Y
    den = one(dtype) - z*Y'
    return num/den
end
lft = linear_fractional_transform



"""
    mobius_transform(z, Y)
    mt(z, Y)

Mobius_transform transforms the upper half complex plane to the unit disk in the complex plane, in which point `Y` is mapped to the center of `D`
"""
function mobius_transform(z, Y) 
    #if abs(Y) == 0 @warn "Y=0.0 maps all points to 0.0!" end
    #if imag(z) < 0 @warn "Im(z) ≥ 0 is required."  end
    #if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required."  end
    num = z - Y
    den = z - Y'
    return num/den
end
mt = mobius_transform


"""
    inverse_mobius_transform(z, Y)
    imt(z, Y)
    
Inverse Mobius transform transforms the unit disk in the complex plane to the upper half complex plane.
"""
function inverse_mobius_transform(z, Y)
    #if abs(z) ≥ 1 @warn "|z|<1 is required." end
    #if abs(Y) == 0 @warn "Y=0.0 maps all points to 0!" end
    #if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required." end
    dtype = eltype(z)
    num = Y - z*Y'
    den = one(dtype) - z
    return num/den
end
imt = inverse_mobius_transform





