using NevanlinnaAC
path = "/home/sliang/JuliaCode/NevanlinnaAC/test/F03"
cd(path)

ngrid = 20
β = 10
otype = Fermi

option = Options(; ngrid, otype)
setprecision(option.precision)

rd = read_to_RawData(path*"/giwn.txt", option)
nd = toNevData(rd, option)

isvalid(nd)
issolvable(nd)

wmesh, Aw = spectral_function(option, rd)