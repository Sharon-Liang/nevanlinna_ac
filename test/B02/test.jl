#ising
using NevanlinnaAC
path = "/home/sliang/JuliaCode/NevanlinnaAC/test/B02"
cd(path)

ngrid = 20
otype = Bose

option = Options(; ngrid, otype)
setprecision(option.precision)


rd = read_to_RawData("/data/sliang/JuliaCMPO/XXZ/Jz_0.00_Jxy_1.00_wid_01/bondD_16_Correlation_Masubara_Freq/giwn_szsz_beta_20.00.txt", option, skipstart=2)
nd = toNevData(rd, option)

isvalid(nd)
issolvable(nd)

wmesh, Aw = spectral_function(option, rd)