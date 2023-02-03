"""
    pick_matrix(x, y)

Construct the Pick matrix of data set `{z, f(z)}` within a unit circle in the complex plane.
"""
function pick_matrix(x::AbstractVector, y::AbstractVector)
    @assert (all(abs.(x) .< 1) && all(abs.(y) .< 1))
    @assert length(x) == length(y) "DimensionMismatch"

    N = length(y)
    pick = similar(y, N, N)
    for i = 1: N, j = 1: N
        num = 1 - y[i] * y[j]'
        den = 1 - x[i] * x[j]'
        pick[i,j] = num/den
    end
    return pick
end


#=
### *Validity*
=#
"""
    isvalid(d::AbstractData)

Determin if the input data ``d`` within its definition domain.
"""
isvalid(d::AbstractData) = "Cannot determin." 


"""
    isvalid(d::GenSchurData)

Determin if ``d.grid`` is in the upper half complex plain and ``d.value`` with in a unit circle centered at 0.0.
"""
function isvalid(d::GenSchurData)
    return all(imag.(d.grid) .≥ 0) &&  all(abs.(d.value) .≤ 1)
end


"""
    isvalid(d::NevData)

Determin if ``d.grid`` and ``d.value`` are all in the upper half complex plain.
"""
function isvalid(d::NevData)
    return all(imag.(d.grid) .≥ 0) &&  all(abs.(d.value) .≥ 0)
end



#=
### *Solvability*
=#

"""
    issolvable(d::AbstractData)

Determin if the input data ``d`` is solvable by corresponding extra interpolation algorithm.
"""
issovable(d::AbstractData) = "Cannot determin." 


"""
    issovable(d::GenSchurData[; tolerance = 1.e-10])
    
Check if the initial data `{z, f(z)}` can be interpolated by generalized schur algorithm.
"""
function issovable(d::GenSchurData; tolerance::AbstractFloat = 1.e-10, show_warning::Bool=false)
    @assert isvalid(d)
        
    x = _mti.(d.grid; show_warning)
    pick = pick_matrix(x, d.value)

    evals = eigvals(pick) 
    return all((real.(evals) .+ tolerance) .>= 0), minimum(real.(evals))
end


"""
    issolvable(d::NevData[; tolerance=1.e-10])
    
Check if the initial data `{z, f(z)}` can be interpolated by generalized schur algorithm for Nevanlinna functions.
"""
function issolvable(d::NevData; tolerance::AbstractFloat = 1.e-10, show_warning::Bool=false)
    @assert isvalid(d)

    x = _mti.(d.grid; show_warning)
    y = _mti.(d.value; show_warning)
        
    pick = pick_matrix(x, y)
    evals = eigvals(pick) 
    return all((real.(evals) .+ tolerance) .>= 0.), minimum(real.(evals))
end


#=
### *Interpolation* : *help functions*
=#
"""
    _recursion(A::AbstractMatrix, θn::Number) 

Recursion relation of contractive functions `θ(n)` and `θ(n-1)`,the relation is:
```math
θ(n-1) = [A[1,1]*θ(n) + A[1,2]]/[A[2,1]*θ(n)+A[2,2]]
'''

See also : [`_inv_recursion`])(@ref)
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

See also : [`_recursion`])(@ref)
"""
function _inv_recursion(A::AbstractMatrix, θp::Number)
    num = -A[2,2] * θp + A[1,2]
    den = A[2,1] * θp - A[1,1]
    return num/den
end


"""
    coefficient(z, xj, ϕj) 

Calculate the coefficient of the recursion relation in generalized_schur algorithm.
"""
function coefficient(z::Number, xj::Number, ϕj::Number, show_warning::Bool=false) 
    A = zeros(typeof(xj),2,2)
    A[1,1] = mt(z, xj; show_warning)
    A[1,2] = ϕj
    A[2,1] = ϕj' * mt(z, xj)
    A[2,2] = one(typeof(xj))
    return A
end


#=
### *Interpolation* : *core*
=#

"""
    schur_parameter(d::GenSchurData{T}) where T
    
Evaluate Schur parameters for contractive functions `y(x)` within a unit circle, return a list of Schur parameters.
"""
function schur_parameter(d::GenSchurData, show_warning::Bool=false) 
    x = d.grid
    y = d.value
    T = eltype(x)

    L = lastindex(x)

    # sparam: a verctor to store Schur parameters
    sparam = zeros(T, L)
    sparam[1] = x[1]

    #Mout: a matrix of 2×2×(L-1)×L arrays to store the coefficient matrices in each step
    Mout = zeros(T, 2, 2, L-1, L)

    #Mat: a vector of 2×2 arrays to store continued products of coefficient matrix in each step
    Macc = fill(eye(T,2), L)

    #factor: a vector of 2×2×(L-1) arrays to store the final coefficient matrix in each step
    factor = zeros(T, 2, 2, L-1)

    for j = 1 : (L-1)
        for k = j : L
            prod = coefficient(x[k], x[j], sparam[j], show_warning)
            Macc[k] *= prod
            Mout[:, :, j, k] = prod
        end
        sparam[j+1] = _inv_recursion(Macc[j+1], y[j+1])
        factor[:, :, j] = Macc[j+1]
    end

    @save "./schur_parameter.jld" sparam factor Mout
    return sparam
end


"""
    generalized_schur(z::Number, x::AbstractVector{T}, y::AbstractVector{T}[; init_func::Function = x->zero(T)]) where T

The generalized Schur algorithm that interpolates beween `{x,y}` and generate a contractive function `f(z)`, return its value at `z`.
"""
function generalized_schur(zmesh::AbstractVector{T}, d::GenSchurData{T}; show_warning::Bool=false) where T
    @assert isvalid(d) "Invalid data"

    sparam = schur_parameter(d, show_warning)

    #Cmat: coefficient matrix
    func = z -> begin
        Cmat = eye(T, 2)
        for j in eachindex(d.grid)
            Cmat *= coefficient(z, d.grid[j], sparam[j], show_warning)
        end
        return _recursion(Cmat, zero(T))
    end

    return map(func, zmesh)
end


function generalized_schur(zmesh::AbstractVector{T}, d::GenSchurData{T}, params::AbstractArray{T}; show_warning::Bool=false) where T
    @assert isvalid(d) "Invalid data"

    sparam = schur_parameter(d, show_warning)

    #Cmat: coefficient matrix
    func = z -> begin
        Cmat = eye(T, 2)
        for j = 1 : length(d.grid)
            Cmat *= coefficient(z, x[j], sparam[j], show_warning)
        end
        return _recursion(Cmat, hardy_expand(z, params))
    end
    return map(func, zmesh)
end


"""
    nevanlinna(z::T, d::NevData{T}, args...; show_warning::Bool = false) where T

The Nevanlinna Interpolation algorithm that interpolates beween `{x,y}` and generate a nevanlinna function `f(z)`, return its value at `z`.
    
"""
function nevanlinna(zmesh::AbstractVector{T}, d::NevData{T}, args...; show_warning::Bool = false) where T
    @assert isvalid(d) "Invalid data"

    res = generalized_schur(zmesh, toGenSchurData(d), args...; show_warning)

    return _imti.(res)
end


"""
    spectral_function(option::Options, d::RawData, args...)

Calculate spectral functions
"""
function spectral_function(option::Options, d::RawData, args...)
    @unpack otype, η = option

    nd = toNevData(d, option)
    mesh = make_mesh(option)
    wmesh = real.(mesh)

    #Construt the Nevanlinna function
    Nz = nevanlinna(mesh, nd, args...)

    #convert to spectral function A(ω) = -ImGz/π
    if otype == Fermi
        Aw = imag(Nz)/π 
    else
        num = @. wmesh * imag(Nz) - η * real(Nz)
        den = @. (wmesh^2 + η^2)*π
        Aw = @. num/den
    end

    open("./spectral_function.txt", "w") do io
        write(io, "              ω                     A(ω)             \n")
        write(io, "------------------------     ------------------------\n")
        writedlm(io, [wmesh Aw])
    end
    return wmesh, Aw   
end