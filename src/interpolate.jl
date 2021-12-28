
"""
    Conformal transforms used
"""
function mti(z::Ctype) 
    dtype = typeof(z)
    return mt(z, one(dtype)im)
end

function imti(z::Ctype)
    dtype = typeof(z)
    return imt(z, one(dtype)im)
end


"""
    Pick matrix of initial data {z, f(z)} with in a unit cell
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
        x = mti.(x)
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
        x = mti.(x)
        y = mti.(y)
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
function coefficient(z::T, xj::T, ϕj::T) where T<:Ctype
    dtype = typeof(z)
    A = zeros(dtype,2,2)
    A[1,1] = mt(z, xj)
    A[1,2] = ϕj
    A[2,1] = ϕj'* mt(z, xj)
    A[2,2] = one(dtype)
    return A
end


"""
    core: evaluate 'Schur parameters' for contractive functions
    y1 within a unit circle
"""
function schur_parameter(x::AbstractVector, y::AbstractVector)
    M = length(y) ; dtype = eltype(x)
    ϕ = zeros(dtype, M); ϕ[1] = y[1]
    abcd = [eye(dtype,2) for i=1:M]
    for j = 1:(M-1)
        for k=j:M
            prod = coefficient(x[k], x[j], ϕ[j])
            abcd[k] *= prod
        end
        ϕ[j+1] = inv_recursion(abcd[j+1], y[j+1])
        #if abs(ϕ[j+1]) ≥ 1.0 
        #    msg = @sprintf "%i-th Schur parameter ≥ 1 with absolute value: %.5f" j+1 abs(ϕ[j+1])
        #    @warn msg
        #end
    end
    return ϕ
end

"""
    Generalized Schur algorithm
"""
function generalized_schur(z::Number, x::AbstractVector, y::AbstractVector;
    optim::Symbol = :none)
    if all(imag.(x) .≥ 0) == false 
        @warn "Initial data should be in the upper half complex plane"
    elseif all(abs.(y) .≤ 1) == false 
        @error "Target data should be in the unit circle"
    else
        M = length(y); dtype = eltype(y)
        z = dtype(z)
        ϕ = schur_parameter(x,y)
        abcd = eye(dtype,2)
        for j = 1:M
            abcd *= coefficient(z,x[j],ϕ[j])
        end

        if optim == :none 
            θm(z::Number) = zero(dtype)
        end
        return recursion(abcd, θm(z))
    end
end

function generalized_schur(z::AbstractArray, x::AbstractVector, y::AbstractVector;
    optim::Symbol = :none)
    res = similar(z, Ctype)
    for i=1:length(z)
        res[i] = generalized_schur(z[i], x, y, optim=optim)
    end
    return res
end

"""
    Nevanlinna Interpolation algorithm
"""
function nevanlinna(z::Number, x::AbstractVector, y::AbstractVector;
    optim::Symbol = :none)
    if all(imag.(x) .≥ 0) == false 
        @warn "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) .≥ 0) == false
        @warn "Target data should be in the upper half complex plane"
    end
    y = mti.(y)
    res = generalized_schur(z, x, y, optim=optim)
    return imti(res)
end

function nevanlinna(z::AbstractArray, x::AbstractVector, y::AbstractVector;
    optim::Symbol = :none)
    res = similar(z, Ctype)
    for i=1:length(z)
        res[i] = nevanlinna(z[i], x, y, optim=optim)
    end
    return res
end