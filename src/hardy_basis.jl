"""
    H2basis(z, k)

The `k`-th orthonormal basis of `H2` space.
"""
function H2basis(z::Number, k::Int64)
    return 1/(√π*(z+oneunit(z)im)) * _mti(z)^k
end


"""
    hardy_expand(z, Nh::Int64, params::AbstractArray)

    generate a function `f(z)` by the Hardy basis up to order ``div(length(params), 2)``, the corresponding coefficients are stored in array `params`. 
"""
function hardy_expand(z::Number, params::AbstractArray)
    N = div(length(params), 2)
    res = [params[2n-1]*H2basis(z,n-1) + params[2n]*conj(H2basis(z,n-1)) for n=1:N] |> sum
    return res
end


