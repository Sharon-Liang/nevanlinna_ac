using NevanlinnaAC, DelimitedFiles

path = "/home/sliang/JuliaCode/NevanlinnaAC/test/F13"
cd(path)

ngrid = 20
otype = Fermi

option = Options(; ngrid, otype, precision=128)
setprecision(option.precision)

rd = read_to_RawData(path*"/giwn.txt", option)
nd = toNevData(rd, option)

isvalid(nd)
issolvable(nd)


wmesh, Aw = spectral_function(option, rd)
open("./spectral_function_init.txt", "w") do io
    write(io, "              ω                     A(ω)             \n")
    write(io, "------------------------     ------------------------\n")
    writedlm(io, [wmesh Aw])
end

using Optim, Zygote
nhardy = 10
p0 = zeros(2*nhardy)
l = ps -> loss(ps, rd, option, λ = 1.e-4)

#Zygote.gradient(l, p0)[1]

res = Optim.optimize(l, gradient_function(l, p0), p0, LBFGS(), Optim.Options(show_trace=true, show_every=10))
wmesh_opt, Aw_opt = spectral_function(option, rd, res.minimizer)
open(path*"/spectral_function_opt.txt", "w") do io
    write(io, "              ω                     A(ω)             \n")
    write(io, "------------------------     ------------------------\n")
    writedlm(io, [wmesh_opt Aw_opt])
end