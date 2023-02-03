#ising
using NevanlinnaAC
path = "/home/sliang/JuliaCode/NevanlinnaAC/test/B01"
cd(path)

ngrid = 20
β = 10
otype = Bose

option = Options(; ngrid, otype)
setprecision(option.precision)


rd = read_to_RawData("/data/sliang/JuliaCMPO/TFIsing/J_1.00_G_1.00_wid_01/bondD_16_Correlation_Masubara_Freq/giwn_pzpz_beta_10.00.txt", option, skipstart=2)
nd = toNevData(rd, option)

isvalid(nd)
issolvable(nd)

wmesh, Aw = spectral_function(option, rd)