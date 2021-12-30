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
    hcpp = tocomplex(hiwn_cpp[:,2], hiwn_cpp[:,3])

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
    
    xcpp = one(eltype(hiwn_cpp))im * hiwn_cpp[:,1]


    #compare spectrum
    A0 = readdlm("./data/gaussian/A_gaussian_1.txt")
    A1 = readdlm("./data/gaussian/cpp/spectrum/A_gaussian_1_N_21_F64.txt")
    A2 = readdlm("./data/gaussian/cpp/spectrum/A_gaussian_1_N_21_rand_1_F64.txt")
    A3 = readdlm("./data/gaussian/cpp/spectrum/A_gaussian_1_N_21_rand_2_F64.txt")

    A11 = readdlm("./data/gaussian/cpp/spectrum/A_gaussian_1_N_21_F64.txt")
    A21 = readdlm("./data/gaussian/cpp/spectrum/A_gaussian_1_N_21_rand_1_F64.txt")
    A31 = readdlm("./data/gaussian/cpp/spectrum/A_gaussian_1_N_21_rand_2_F64.txt")

    plot(A0[:,1], A0[:,2], label="theory", line=(:black,2))
    plot!(A1[:,1], A1[:,2],label="c++", line=(:dash,2))
    plot!(A2[:,1], A2[:,2],label="c++ rand 1", line=(:dash,2))
    plot!(A3[:,1], A3[:,2],label="c++ rand 2", line=(:dash,2))

    plot!(A1[:,1], A1[:,2],label="c++", line=(:dot,2))
    plot!(A2[:,1], A2[:,2],label="c++ rand 1", line=(:dot,2))
    plot!(A3[:,1], A3[:,2],label="c++ rand 2", line=(:dot,2))

    plot!(xlim=(-5,5))
    savefig("./compare_A_F128.pdf")

    h1 = d1[:,2] .+ d1[:,3] .* 1.0im
    hj = mt.(-y1, one(ComplexF64))


    setprecision(BigFloat, 1024)
    plot()
    for ftype in [Float64]
        tname = @sprintf "%s" ftype
        s, abcd, theta = schur_parameter(ftype, x1n, y1n)
        sc, abcdc, thetac = schur_parameter(ftype, xcpp, hcpp)
        plot!(imag.(xcpp), abs.(s .- sc), line=(:dash,1), 
            marker=(:circle,5,stroke(0.3)), 
            label = tname)
    end
    plot!(xlable="ωn", ylabel="|s(x:julia) - sc(x:c++)|", yaxis = :log)
    savefig("./schur_parameter_loss.pdf")

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



