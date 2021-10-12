struct Giωn
    ωn::Real
    val::Complex
end

struct MasubaraGF
    number::Int64
    GF::Vector{Giωn}
end

function make_input(path::String, num::Int64)
    data = readdlm(path)
    fnum = min(num, length(data[:,1])) # number of frenquencies used
    res = [Giωn(data[i,1],data[i,2]+1.0im*data[i,3]) for i = 1:fnum]
    return MasubaraGF(fnum, res |> reverse)
end

function make_input(path::String)
    data = readdlm(path)
    fnum = length(data[:,1]) # number of frenquencies used
    res = [Giωn(data[i,1],data[i,2]+1.0im*data[i,3]) for i = 1:fnum]
    return MasubaraGF(fnum, res |> reverse)
end

function pick_matrix(G::MasubaraGF)
    # input :freq: Masubara frenquencies; val: G
    N = G.number; g = G.GF
    pick = zeros(ComplexF64,N, N)
    for i = 1:N, j = 1:N
        up = 1 - λ(g[i])*λ(g[j])'
        down = 1 - h(1.0im*g[i].ωn)*h(1.0im*g[j].ωn)'
        pick[i,j] = up/down
    end
    return pick
end

function isNevanlinnasolvable(G::MasubaraGF; err::Float64 = 1.e-6)
    # check if the pick_matrix is semi positive semidefinite
    pick = pick_matrix(G)
    evals = eigvals(pick) 
    return all((evals .+ err) .>= 0), minimum(evals)
end

function λ(G::Giωn)
    return h(-G.val)
end

function θk(fac::Matrix, θp::Number)
    # θp = θ(k-1)
    up = -fac[2,2] * θp + fac[1,2]
    down = fac[2,1] * θp - fac[1,1]
    return up/down
end

function core(G::MasubaraGF; etype::DataType=Float64)
    dtype = Complex{etype}
    M = G.number; g = G.GF
    phis = zeros(dtype, M) #store λk(iω_(k+1))
    phis[1] = λ(g[1])
    abcds = Vector{Matrix{dtype}}(undef,M)
    for i = 1:M abcds[i] = eye(dtype,2) end
    for j = 1:(M-1)
        for k=j:M
            prod = zeros(dtype,2,2)
            prod[1,1] = h1(1.0im*g[k].ωn,1.0im*g[j].ωn)
            prod[1,2] = phis[j]
            prod[2,1] = phis[j]' * h1(1.0im*g[k].ωn,1.0im*g[j].ωn)
            prod[2,2] = 1. + 0.0im
            abcds[k] *= prod
        end
        phis[j+1] = θk(abcds[j+1], λ(g[j+1]))
    end
    return phis
end

function evaluation(z::Number, G::MasubaraGF; 
    optim::Symbol = :none, etype::DataType=Float64)
    dtype = Complex{etype}
    M = G.number
    g = G.GF
    phis = core(G, etype=etype)
    abcd = eye(dtype, 2)
    for j = 1:M
        prod = zeros(dtype, 2,2)
        prod[1,1] = h1(z,1.0im*g[j].ωn)
        prod[1,2]= phis[j]
        prod[2,1] = phis[j]' * h1(z,1.0im*g[j].ωn)
        prod[2,2] = 1. + 0.0im
        abcd *= prod
    end
    
    if optim == :none
        θm(z::Number) = 0.0 + 0.0im
    end
    
    up = abcd[1,1] * θm(z) + abcd[1,2]
    down = abcd[2,1] * θm(z) + abcd[2,2]
    θ = up/down
    return invh(θ)
end

function spectrum_density(ω::Real, η::Real, G::MasubaraGF;
    etype::DataType=Float64)
    z = ω + 1.0im * η
    return 2*(evaluation(z, G, etype=etype) |> imag)
end

