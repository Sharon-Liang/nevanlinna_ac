using nevanlinna_ac
using DelimitedFiles
using Printf
using DoubleFloats
using Debugger
using Plots; gr(xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
default(palette = palette(:okabe_ito))

function tomatrix(mv::AbstractVector)
    M = length(mv)
    T = eltype(real.(mv[1])) 
    res = zeros(T, M, 8)
    for i = 1:M, j = 1:4
        A = mv[i] |> transpose
        res[i,2j-1] = A[j] |> real
        res[i,2j] = A[j] |> imag
    end
    return res
end

function tocomplex(a::Vector{T}, b::Vector{T}) where T
    I = one(T)im
    res = a .+ I .* b
    return res
end

# n in [5, 21, 40]
n = 21
prec = 64
#function ploss(preci::Number, n::Integer)
#    setprecision(preci)
    
    # c++ data
    hiwn_path = @sprintf "./data/gaussian/cpp/F%i/hiwn_gaussian_1_N_%i.txt" prec n
    sparam_path = @sprintf "./data/gaussian/cpp/F%i/sparam_gaussian_1_N_%i.txt" prec n
    abcd_path = @sprintf "./data/gaussian/cpp/F%i/abcd_gaussian_1_N_%i.txt" prec n
    theta_path = Vector{String}(undef, n-1)
    for i = 1:n-1
        theta_path[i] = @sprintf "./data/gaussian/cpp/F%i/theta_%i_gaussian_1_N_%i.txt" prec i n
    end

    hiwn_cpp = readdlm(hiwn_path)
    hcpp = tocomplex(hiwn_data[:,2], hiwn_data[:,3])

    sparam_data = readdlm(sparam_path)
    scpp = tocomplex(sparam_data[:,2], sparam_data[:,3])

    abcd_cpp = readdlm(abcd_path)
    theta_cpp = [readdlm(theta_path[i]) for i=1:n-1]
    
    #julia data
    giwn_path = "./data/gaussian/giwn_gaussian_1.txt"
    x1, y1 = readGF(giwn_path, num=n)
    x1n, y1n = toGeneralizedSchurdata(x1, y1, :f)
    #confirm x1n = hiwn_data[:,1]
    #maximum(abs.(y1n .- hc)) 2.220462989844635e-16 = eps(Float64)
    
    s, abcd, theta = schur_parameter(x1n, y1n)

    abcd = tomatrix(abcd)
    theta = [tomatrix(theta[i,:]) for i=1:n-1]


    ds = s .- scpp

    dθ = similar(theta_cpp)
    for i = 1:n-1
        dθ[i] = similar(theta_cpp[i])
        dθ[i][:,1] = theta_cpp[i][:,1]
        dθ[i][:,2:end] = abs.(theta[i][i:end,:].- theta_cpp[i][:,2:end])
    end

    max_da = [maximum(abs.(dθ[i][:,2:3])) for i=1:n-1]
    max_db = [maximum(abs.(dθ[i][:,4:5])) for i=1:n-1]
    max_dc = [maximum(abs.(dθ[i][:,6:7])) for i=1:n-1]
    max_dd = [maximum(abs.(dθ[i][:,8:9])) for i=1:n-1]
#    return sp[:,1], abs.(dimg), abs.(fac_out .- abcd)
#end


function plot_fac_diff(x1::AbstractVector, fac::AbstractMatrix)
    x1 = x1[1:end-1]
    plot(x1, fac[:,2], marker=(:circle, 3), line=(:dash,1), label="a.real")
    plot!(x1, fac[:,3], marker=(:circle, 3), line=(:dash,1), label="a.imag")
    plot!(x1, fac[:,4], marker=(:circle, 3), line=(:dash,1), label="b.real")
    plot!(x1, fac[:,4], marker=(:circle, 3), line=(:dash,1), label="a.real")
end



