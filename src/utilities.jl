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
function linear_fractional_transform(z::Number, Y::Number)
    if abs(z) ≥ 1 || abs(Y) ≥ 1
        @error "|z|<1 and |Y|<1 are required."
    else
        num = z - Y
        den = 1 - z*Y'
        return num/den
    end
end
lft = linear_fractional_transform

"""
    Mobius_transform transforms the upper half complex plane
    to the unit disk in the complex plane
"""
function mobius_transform(z::Number, Y::Number)
    if abs(Y) == 0
        @error "Y should not be 0!"
    elseif imag(z) < 0 || imag(Y) < 0
        @error "Im(z) ≥ 0 and Im(Y) ≥ 0 are required."  
    else
        num = z - Y
        den = z - Y'
        num == den ? (return 1) : (return num/den)
    end
end
mt = mobius_transform


"""
    Inverse Mobius transform transforms the unit disk in the complex plane
    to the upper half complex plane.
"""
function inverse_mobius_transform(z::Number, Y::Number)
    if abs(z) ≥ 1
        @error "|z|<1 is required."
    elseif abs(Y) == 0 
        @error "Y should not be 0!"
    elseif imag(Y) < 0
        @error "Im(Y) ≥ 0 is required."  
    else
        num = Y - z*Y'
        den = 1 - z
        return num/den
    end
end
imt = inverse_mobius_transform





