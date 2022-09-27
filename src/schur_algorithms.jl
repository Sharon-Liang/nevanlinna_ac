"""
    _mti(z)

mobius_transform(z, 1im)
"""
_mti(z) = mt(z, oneunit(z)im)

"""
    _imti(z)

inverse_mobius_transform(z, 1im)
"""
_imti(z) = imt(z, oneunit(z)im)


"""
    ispossemidef(A::AbstractMatrix)

Check if a matrix `A` is positive semidefinite.
"""
function ispossemidef(A::AbstractMatrix)
    evals = eigvals(A)
    return all(evals .>= 0)
end


"""
    pick_matrix(x, y)

Pick matrix of initial data `{z, f(z)}` with in a unit cell.
"""
function pick_matrix(x::AbstractVector, y::AbstractVector)
    if (all(abs.(x) .> 0) && all(abs.(y) .> 0)) == false
        @error "pick_matrix: initial data and target data should be in a unit cell"
    elseif length(x) != length(y)
        @error DimensionMismatch
    else
        N = length(y)
        pick = similar(y, N, N)
        for i = 1: N, j = 1: N
            num = 1 - y[i] * y[j]'
            den = 1 - x[i] * x[j]'
            pick[i,j] = num/den
        end
        return pick
    end
end


"""
    isGeneralizedSchursovable(x, y[;tolerance=1.e-10])
    
Check if the initial data `{z, f(z)}` can be interpolated by generalized schur algorithm.
"""
function isGeneralizedSchursovable(x::AbstractVector, y::AbstractVector; 
    tolerance::AbstractFloat = 1.e-10)
    if all(imag.(x) .≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(abs.(y) .≤ 1) == false 
        @error "Target data should be in the unit circle"
    else
        x = _mti.(x)
        pick = pick_matrix(x, y)
        evals = eigvals(pick) 
        return all((real.(evals) .+ tolerance) .>= 0), minimum(real.(evals))
    end
end


"""
    isNevanlinnasolvable(x, y[;tolerance=1.e-10])
    
Check if the initial data `{z, f(z)}` can be interpolated by generalized schur algorithm for Nevanlinna functions.
"""
function isNevanlinnasolvable(x::AbstractVector, y::AbstractVector; 
    tolerance::AbstractFloat = 1.e-10)
    if all(imag.(x) .≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) .≥ 0) == false
        @error "Target data should be in the upper half complex plane"
    else
        x = _mti.(x)
        y = _mti.(y)
        pick = pick_matrix(x, y)
        evals = eigvals(pick) 
        return all((real.(evals) .+ tolerance) .>= 0.), minimum(real.(evals))
    end
end


"""
    _recursion(A::AbstractMatrix, θn::Number) 

Recursion relation of contractive functions `θ(n)` and `θ(n-1)`,the relation is:
```math
θ(n-1) = [A[1,1]*θ(n) + A[1,2]]/[A[2,1]*θ(n)+A[2,2]]
'''
"""
function _recursion(A::AbstractMatrix, θn::Number)
    num = A[1,1] * θn + A[1,2]
    den = A[2,1] * θn + A[2,2]
    return num/den
end


"""
    _inv_recursion(A::AbstractMatrix, θn::Number) 

Recursion relation of contractive functions `θ(n)` and `θ(n+1)`,the relation is:
```math
θ(n+1) = [-A[2,2]*θ(n) + A[1,2]]/[A[2,1]*θ(n)-A[1,1]]
'''
"""
function _inv_recursion(A::AbstractMatrix, θp::Number)
    num = -A[2,2] * θp + A[1,2]
    den = A[2,1] * θp - A[1,1]
    return num/den
end


"""
    _coefficient(z, xj, ϕj) 

Calculate the coefficient of the recursion relation in generalized_schur algorithm.
"""
function coefficient(z::Number, xj::Number, ϕj::Number) 
    A = zeros(typeof(z),2,2)
    A[1,1] = mt(z, xj)
    A[1,2] = ϕj
    A[2,1] = ϕj' * mt(z, xj)
    A[2,2] = one(typeof(z))
    return A
end


"""
    schur_parameter([ftype::DataType = T,] x::AbstractVector{T}, y::AbstractVector{T}) where T
    
Evaluate Schur parameters for contractive functions `y(x)` within a unit circle, return a list of Schur parameters
"""
function schur_parameter(x::AbstractVector{T}, y::AbstractVector{T}) where T
    M = length(y)
    ϕ = zeros(T, M); ϕ[1] = y[1]
    abcd = fill(eye(T,2), M)
    abcd_out = fill(zeros(T,2,2), M-1, M)
    factor = fill(zeros(T,2,2), M-1)
    for j = 1:(M-1)
        for k=j:M
            prod = coefficient(x[k], x[j], ϕ[j])
            abcd[k] *= prod
            abcd_out[j,k] = prod
        end
        ϕ[j+1] = _inv_recursion(abcd[j+1], y[j+1])
        factor[j] = abcd[j+1]
        #if abs(ϕ[j+1]) ≥ 1.0 
        #    msg = @sprintf "%i-th Schur parameter ≥ 1 with absolute value: %.5f" j+1 abs(ϕ[j+1])
        #    @warn msg
        #end
    end
    return ϕ
    #return ϕ, factor, abcd_out
end


"""
    generalized_schur(z::Number, x::AbstractVector{T}, y::AbstractVector{T}[; init_func::Function = x->zero(T)]) where T

The generalized Schur algorithm that extrapolates beween `{x,y}` and generate a contractive function `f(z)`, return its value at `z`.
"""
function generalized_schur(z::Number, x::AbstractVector, y::AbstractVector; init_func::Function = z -> zero(eltype(y))) 
    if all(imag.(x) .≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(abs.(y) .≤ 1) == false 
        @error "Target data should be in the unit circle"
    else
        M = length(y)
        ϕ = schur_parameter(x,y)
        abcd = eye(eltype(y),2)
        for j = 1:M
            abcd *= coefficient(z,x[j],ϕ[j])
        end
        return _recursion(abcd, init_func(z))
    end
end


"""
    nevanlinna(z::Number, x::AbstractVector{T}, y::AbstractVector{T}[; init_func::Function = x->zero(T)]) where T

The Nevanlinna Interpolation algorithm that extrapolates beween `{x,y}` and generate a nevanlinna function `f(z)`, return its value at `z`.
    
"""
function nevanlinna(z::Number, x::AbstractVector{T}, y::AbstractVector{T}; init_func::Function = z -> zero(T)) where T
    if all(imag.(x) .≥ 0) == false 
        @warn "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) .≥ 0) == false
        @warn "Target data should be in the upper half complex plane"
    end
    y = _mti.(y)
    res = generalized_schur(z, x, y; init_func)
    return _imti(res)
end
