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
function ispossemidef(A::Matrix{T} where T<:Number)
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
    Check if the Pick matrix is positive semidefinite
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
        return all((evals .+ tolerance) .>= 0), minimum(evals)
    end
end

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
        return all((evals .+ tolerance) .>= 0), minimum(evals)
    end
end


"""
    Recursion relation of contractive functions θp and θn
    θp : previous θ(z)
    θn : next θ(z)
    recursion: θp = (a*θn + b)/(c*θn+d)
    inv_recursion: θn = (-d*θp + b)/(c*θp - a)
    The target function is θ_0 
"""
function recursion(A::Matrix{T}, θn::T) where T
    num = A[1,1] * θn + A[1,2]
    den = A[2,1] * θn + A[2,2]
    return num/den
end

function inv_recursion(A::Matrix{T}, θp::T) where T
    num = -A[2,2] * θp + A[1,2]
    den = A[2,1] * θp - A[1,1]
    return num/den
end


"""
    coefficient of the recursion relation in generalized_schur algorithm
"""
function coefficient(z::T, xj::T, ϕj::T) where T
    A = zeros(T,2,2)
    A[1,1] = mt(z, xj)
    A[1,2] = ϕj
    A[2,1] = ϕj' * mt(z, xj)
    A[2,2] = one(T)
    return A
end


"""
    core: evaluate 'Schur parameters' for contractive functions
    y1 within a unit circle
    outpus: ϕ: Schur_parameters
    factor[i]: abcd of i-the contractive function
    abcd_out[i,j]: abcd of i-th contractive function at point x[j]
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
        ϕ[j+1] = inv_recursion(abcd[j+1], y[j+1])
        factor[j] = abcd[j+1]
        #if abs(ϕ[j+1]) ≥ 1.0 
        #    msg = @sprintf "%i-th Schur parameter ≥ 1 with absolute value: %.5f" j+1 abs(ϕ[j+1])
        #    @warn msg
        #end
    end
    return ϕ
    #return ϕ, factor, abcd_out
end


function schur_parameter(ftype::DataType, x::AbstractVector{T}, y::AbstractVector{T}) where T
    ctype = Complex{ftype}
    x = ctype.(x)
    y = ctype.(y)
    return schur_parameter(x,y)
end


"""
    Generalized Schur algorithm
"""
function generalized_schur(z::T, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    if all(imag.(x) .≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(abs.(y) .≤ 1) == false 
        @error "Target data should be in the unit circle"
    else
        M = length(y)
        ϕ = schur_parameter(x,y)
        abcd = eye(T,2)
        for j = 1:M
            abcd *= coefficient(z,x[j],ϕ[j])
        end
        return recursion(abcd, init_func(z))
    end
end

function generalized_schur(ftype::DataType, z::T, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    ctype = Complex{ftype}
    z = ctype(z)
    x = ctype.(x)
    y = ctype.(y)
    return generalized_schur(z, x, y; init_func)
end

function generalized_schur(z::AbstractArray{T}, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    res = similar(z, T)
    for i=1:length(z)
        res[i] = generalized_schur(z[i], x, y; init_func)
    end
    return res
end

function generalized_schur(ftype::DataType, z::AbstractArray{T}, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    ctype = Complex{ftype}
    z = ctype.(z)
    x = ctype.(x)
    y = ctype.(y)
    return generalized_schur(z, x, y; init_func)
end


"""
    Nevanlinna Interpolation algorithm
"""
function nevanlinna(z::T, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    if all(imag.(x) .≥ 0) == false 
        @warn "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) .≥ 0) == false
        @warn "Target data should be in the upper half complex plane"
    end
    y = _mti.(y)
    res = generalized_schur(z, x, y; init_func)
    return _imti(res)
end

function nevanlinna(ftype::DataType, z::T, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    ctype = Complex{ftype}
    z = ctype(z)
    x = ctype.(x)
    y = ctype.(y)
    return nevanlinna(z, x, y; init_func)
end

function nevanlinna(z::AbstractArray{T}, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    res = similar(z, T)
    for i=1:length(z)
        res[i] = nevanlinna(z[i], x, y; init_func)
    end
    return res
end

function nevanlinna(ftype::DataType, z::AbstractArray{T}, x::AbstractVector{T}, y::AbstractVector{T};
    init_func::Function = z -> zero(T)) where T
    ctype = Complex{ftype}
    z = ctype.(z)
    x = ctype.(x)
    y = ctype.(y)
    return nevanlinna(z, x, y; init_func)
end