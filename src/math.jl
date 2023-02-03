"""
    Identity matrix
"""
function eye(n::Int64)
    return Matrix{Float64}(I,n,n)
end

function eye(dtype::DataType, n::Int64)
    return Matrix{dtype}(I,n,n)
end


#=
### *Conformal transforms*
=#

"""
    linear_fractional_transform(z, Y)
    lft(z, Y)
    
Linear fractional transform transforms the close unit disk `D` in the complex plane to itself, in which point `Y` is mapped to the center of `D`. Moreover, it is a one to one mapping of the boundary of `D`.
"""
function linear_fractional_transform(z, Y; show_warning::Bool=false) 
    if show_warning
        if abs(z) ≥ 1 @warn "|z|<1 is required." end
        if abs(Y) ≥ 1 @warn "New center |Y|<1 is required." end
    end
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
function mobius_transform(z, Y; show_warning::Bool=false) 
    @assert abs(Y) != 0 
    if show_warning
        if  @warn "Y=0.0 maps all points to 0.0!" end
        if imag(z) < 0 @warn "Im(z) ≥ 0 is required."  end
        if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required."  end
    end
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
function inverse_mobius_transform(z, Y; show_warning::Bool=false)
    @assert abs(Y) != 0 
    if show_warning
        if abs(z) ≥ 1 @warn "|z|<1 is required." end
        if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required." end
    end
    dtype = eltype(z)
    num = Y - z*Y'
    den = one(dtype) - z
    return num/den
end
imt = inverse_mobius_transform

"""
    _mti(z)

mobius_transform(z, 1im)
"""
_mti(z; show_warning::Bool=false) = mt(z, oneunit(z)im; show_warning)


"""
    _imti(z)

inverse_mobius_transform(z, 1im)
"""
_imti(z; show_warning::Bool=false) = imt(z, oneunit(z)im; show_warning)


#=
### *Hardy basis H2*
=#
"""
    H2basis(z, k)

The `k`-th orthonormal basis of `H2` space.
"""
function H2basis(z::Number, k::Int64)
    return 1/(√π*(z+oneunit(z)im)) * _mti(z)^k
end


"""
    hardy_expand(z, params::AbstractVector)

    generate a function `f(z)` by the Hardy basis up to order ``div(length(params), 2)``, the corresponding coefficients are stored in vector `params`. 
"""
function hardy_expand(z::Number, params::AbstractArray)
    N = div(length(params), 2)
    res = [params[2n-1]*H2basis(z,n-1) + params[2n]*conj(H2basis(z,n-1)) for n=1:N] |> sum
    return res
end


