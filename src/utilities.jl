function eye(n::Int64)
    res = zeros(n,n)
    for i = 1:n res[i,i] = 1.0 end
    return res
end

function eye(dtype::DataType, n::Int64)
    res = zeros(dtype, n, n)
    I = 1.0 |> dtype
    for i = 1:n res[i,i] = I end
    return res	
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

function Masubara_GF(n::Int64, A::Function, β::Real;
    type::Symbol = :f, Λ::Float64=100., L::Int64=1000000)
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

function h(z::Number)
    #Conformal mapping(Mobius transform) 
    #\bar{C+} -> \bar{D}        
    return (z - 1.0im)/(z + 1.0im)
end

function invh(z::Number)
    return 1.0im*(1 + z)/(1 - z)
end

function h1(z::Number, Y::Number)
    #Conformal mapping(Mobius transform) 
    #C+ -> D
    return (z - Y)/(z - Y')
end

function invh1(z::Complex, Y::Number)
    return (z*Y'-Y)/(z-1)
end

function ispossemidef(A::Matrix)
    evals = eigvals(A)
    return all(evals .>= 0)
end

