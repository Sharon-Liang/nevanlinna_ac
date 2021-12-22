module nevanlinna_ac
__precompile__()

#setprecision(BigFloat, 128)
#BigFloat defalt precision 256

using LinearAlgebra

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
export toNevanlinnadata,
       spectrum

include("utilities.jl")
include("optim.jl")
include("interpolate.jl")
include("GFdata_manipulation.jl")

end