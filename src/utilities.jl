function pauli(symbol::Symbol)
    if symbol==:x return [0. 1.; 1. 0.]
    elseif symbol==:y return [0. -1im; 1im 0.]
    elseif symbol==:z return [1. 0.; 0. -1.]
    elseif symbol==:+ return [0. 1.; 0. 0.]
    elseif symbol==:- return [0. 0.; 1. 0.]
    else
        error("The input should be :x,:y,:z,:+,:-.")
    end
end

function delta(x::Real; c::Real = 0., η::Real = 0.05)
    num = η
    den = ((x-c)^2 + η^2) * π
    return num/den
end

function gaussian(x::Real, fac::AbstractVector)
    # f(x) = 1/(σ√(2π)) exp(-((x-μ)/σ)^2 /2 )
    μ = fac[1]; σ2 = fac[2]
    σ2 <= 0. ? (@error "σ2 should be positive.") : 
    norm = √(2π * σ2)
    m = (x-μ)^2 / σ2
    return 1/norm*exp(-m/2)
end

function multi_gaussian(N::Int64; ms::Vector{Float64}=[-3.,3.], smax::Real=3.)
    fac = rand(N,2)
    fac[:,1] = rand(N) .* (ms[2] - ms[1]) .+ ms[1]
    fac[:,2] = rand(N) .* smax
    function f(x::Real)
        res = 0.
        for i = 1: N
            res += gaussian(x, fac[i,:])
        end
        res/N
    end
    return f, fac
end

function Masubara_freq(n::Int64, β::Real; type::Symbol= :f)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end

function Masubara_GF(n::Int64, A::Function, β::Real, type::Symbol;
    Λ::Float64=100., L::Int64=1000000)
    # G(iωn) = 1/2π ∫dΩ A(Ω)/(iωn - Ω)
    ωn = Masubara_freq(n, β, type=type)
    dL = 2*Λ/L
    Gn = 0.
    for n = 0:L
        ω = -Λ + dL*n
        Gn += A(ω)/(1.0im*ωn-ω) * dL
    end
    return Gn/(2π)
end

function ispossemidef(A::Matrix)
    evals = eigvals(A)
    return all(evals .>= 0.) 
end

function MobiusTransform(z::ComplexF32)
    return (z - 1.0im)/(z + 1.0im)
end

function invMobiusTransform(z::ComplexF32)
    return 1.0im*(1 + z)/(1 - z)
end

function NG(gval::Vector)
    #input: gval: values of masubara GF
end


function pick_matrix(freq::Vector, gval::Vector)
    # input :freq: Masubara frenquencies; val: G
    dim = length(freq)
    if length(val) != dim
        @error DimensionMismatch("dimension of val must match freq")
    end
    pmat = zeros(ComplexF32, dim, dim)
    for i = 1:dim, j = 1:dim
        num = 1 - val[i]*val[j]'
        den = 1 - MobiusTransform(freq[i])*MobiusTransform(freq[j])'
        pmat[i,j] = num/den
    end
    return pmat
end

function isNevanlinnasolvable(freq::Vector, val::Vector)
    # check if the pick_matrix is semi positive semidefinite
    pmat = pick_matrix(freq, val)
    return pmat |> ispossemidef
end

function isNevanlinnaunique(freq::Vector, val::Vector)
    return pick_matrix(freq, val) |> det |> iszero
end