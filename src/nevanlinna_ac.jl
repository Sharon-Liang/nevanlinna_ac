module nevanlinna_ac
__precompile__()

#setprecision(BigFloat, 128)
#BigFloat defalt precision 256

using LinearAlgebra
using Random; Random.seed!()
using DelimitedFiles

export eye, delta, gaussian, multi_gaussian
export Masubara_freq, Masubara_GF
export h, invh, h1, invh1

export Giωn, MasubaraGF, make_input
export pick_matrix, ispossemidef, isNevanlinnasolvable
export core, evaluation, spectrum_density

#input file format: ωn  Re[G(iωn)]  Im[G(iωn)]
include("utilities.jl")
include("interpolate.jl")

end