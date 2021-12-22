"""
    Conformal transforms used
"""
function mti(z::Number) 
    return mt(z, 1.0im)
end

function imti(z::Number)
    return imt(z, 1.0im)
end


"""
    Pick matrix of initial data {z, f(z)} with in a unit cell
"""
function pick_matrix(x::Vector, y::Vector)
    if (all(abs.(x) > 0) && all(abs.(y) > 0)) == false
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
function isGeneralizedSchursovable(x::Vector, y::Vector; 
    tolerance::Float64 = 1.e-6)
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

function isNevanlinnasolvable(x::Vector, y::Vector; 
    tolerance::Float64 = 1.e-6)
    if all(imag.(x) ≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) ≥ 0) == false
        @error "Target data should be in the upper half complex plane"
    end
    x = mti.(x)
    y = mti.(y)
    pick = pick_matrix(x, y)
    evals = eigvals(pick) 
    return all((evals .+ tolerance) .>= 0), minimum(evals)
end


"""
    Recursion relation of contractive functions θp and θn
    θp : previous θ(z)
    θn : next θ(z)
    θp = (a*θn + b)/(c*θn+d)
    The target function is θ_0 
"""
function recursion(A::Matrix, θn::Number)
    num = A[1,1] * θn + A[1,2]
    den = A[2,1] * θn + A[2,2]
    return num/den
end

function inv_recursion(A::Matrix, θp::Number)
    num = -A[2,2] * θp + A[1,2]
    den = A[2,1] * θp - A[1,1]
end


"""
    coefficient of the recursion relation in generalized_schur algorithm
"""
function coefficient(z::Number, xj::Number, ϕj::Number)
    A = zeros(ComplexF64,2,2)
    A[1,1] = mt(z, xj)
    A[1,2] = ϕj
    A[2,1] = ϕj'* mt(z,xj)
    A[2,2] = 1.0 + 0.0im
    return A
end


"""
    core: evaluate 'Schur parameters' for contractive functions
    {y} within a unit circle
"""
function schur_parameter(x::Vector, y::Vector)
    M = length(y); T = ComplexF64
    ϕ = zeros(T, M); ϕ[1] = y[1]
    abcd = [eye(T,2) for i=1:M]
    for j = 1:(M-1)
        for k=j:M
            prod = coefficient(x[k], x[j], ϕ[j])
            abcd[k] *= prod
        end
        ϕ[j+1] = inv_recursion(abcd[j+1], y[j+1])
    end
    return ϕ
end

"""
    Generalized Schur algorithm
"""
function generalized_schur(z::Number, x::Vector, y::Vector;
    optim::Symbol = :none)
    if all(imag.(x) .≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(abs.(y) .≤ 1) == false 
        @error "Target data should be in the unit circle"
    else
        ϕ = schur_parameter(x,y)
        abcd = eye(ComplexF64,2)
        for j = 1:M
            abcd *= coefficient(z,x[j],ϕ[j])
        end

        if optim == :none 
            θm(z::Number) = 0.0 + 0.0im
        end
        return recursion(abcd, θm(z))
    end
end


"""
    Nevanlinna Interpolation algorithm
"""
function nevanlinna(z::Number, x::Vector, y::Vector;
    optim::Symbol = :none)
    if all(imag.(x) ≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) ≥ 0) == false
        @error "Target data should be in the upper half complex plane"
    end
    y = mti.(y)
    res = generalized_schur(z, x, y, optim=optim)
    return imti(res)
end


