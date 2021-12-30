using LinearAlgebra
using DelimitedFiles
using Printf
using DoubleFloats
using Test

"""
    Identity matrix
"""
function eye(dtype::DataType, n::Int64)
    return Matrix{dtype}(I,n,n)
end


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
Mobius_transform: maps 1im to 0
"""
function mti(ftype::DataType, z::T) where T
    I = one(ftype)im
    z = Complex{ftype}(z)
    return mt(z, I)
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
    median result: factor[i]: abcd of i-the contractive function
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


"""Run:
    x ∈ upper half complex plane (iωn)
    y ∈ unit circle in the complex plane
"""
#load input data: Data type : Float64
#c++ result with the same input: 
#    path = "./scpp_F64.txt", "./scpp_F128.txt"
d = readdlm("./input.txt")
I = one(eltype(d))im
x = I * d[:,1]
y = d[:,2] .+ I .* d[:,3]


# float64
ϕ = schur_parameter(x, y)


# Double64
ϕ1 = schur_parameter(Double64, x, y)

# compare with c++
dc64 = readdlm("./scpp_F64.txt")
ϕc = dc64[:,2] .+ one(eltype(dc64))im .* dc64[:,3]
for i=1:length(ϕ)
    @test isapprox(ϕc[i], ϕ[i], rtol=sqrt(eps(eltype(dc64))))
end
