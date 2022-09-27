"""
    OperatorType Bose Fermi

    Operator type of physical observables
"""
@enum OperatorType Bose Fermi

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
    Check if a matrix if positive semidefinite.
"""
function ispossemidef(A::Matrix{T} where T<:Number)
    evals = eigvals(A)
    return all(evals .>= 0)
end


"""
    Linear fractional transform transforms the unit disk in the 
    complex plane to itself
"""
function linear_fractional_transform(z::T, Y::T) where T<:Number
    if abs(z) ≥ 1 @warn "|z|<1 is required." end
    if abs(Y) ≥ 1 @warn "New center |Y|<1 is required." end
    dtype = eltype(z)
    num = z - Y
    den = one(dtype) - z*Y'
    return num/den
end
lft = linear_fractional_transform


"""
    Mobius_transform transforms the upper half complex plane
    to the unit disk in the complex plane
"""
function mobius_transform(z::T, Y::T) where T<:Number
    if abs(Y) == 0 @warn "Y=0.0 maps all points to 0.0!" end
    if imag(z) < 0 @warn "Im(z) ≥ 0 is required."  end
    if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required."  end
    num = z - Y
    den = z - Y'
    return num/den
end
mt = mobius_transform


"""
    Inverse Mobius transform transforms the unit disk in the complex plane
    to the upper half complex plane.
"""
function inverse_mobius_transform(z::T, Y::T) where T<: Number
    if abs(z) ≥ 1 @warn "|z|<1 is required." end
    if abs(Y) == 0 @warn "Y=0.0 maps all points to 0!" end
    if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required." end
    dtype = eltype(z)
    num = Y - z*Y'
    den = one(dtype) - z
    return num/den
end
imt = inverse_mobius_transform





