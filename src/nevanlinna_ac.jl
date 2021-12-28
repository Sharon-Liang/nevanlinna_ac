module nevanlinna_ac
__precompile__()

using LinearAlgebra
using DelimitedFiles
using Printf
using DoubleFloats

export eye, ispossemidef
export linear_fractional_transform, lft,
       mobius_transform, mt, 
       inverse_mobius_transform, imt
export pick_matrix, 
       isGeneralizedSchursovable,
       isNevanlinnasolvable,
       schur_parameter,
       generalized_schur,
       nevanlinna 
export coefficient,
       recursion, inv_recursion
export readGF, toNevanlinnadata,
       spectrum
export Ftype, Ctype

#setprecision(BigFloat, 128)
#BigFloat defalt precision 256

#const Ftype = BigFloat
const Ftype = Double64
const Ctype = Complex{Ftype}

include("utilities.jl")
include("optim.jl")
include("interpolate.jl")
include("GFdata_manipulation.jl")

end