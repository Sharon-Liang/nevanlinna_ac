module NevanlinnaAC
__precompile__()

#=
### *Using Standard Libraries*
=#

using Parameters
using LinearAlgebra, DoubleFloats
using JLD, DelimitedFiles, Printf
using Zygote#, Optim
using FFTW, HCubature

import Base.isvalid
#=
### *Includes And Exports* : *global.jl*
=#

#=
*Summary* :

Define some type aliases and string constants for the NevanlinnaAC toolkit.

*Members* :

```text
Ftype -> Global used numerical type for float numbers
Ctype -> Global used numerical type for complex numbers

#
__LIBNAME__ -> Name of this julia toolkit.
__VERSION__ -> Version of this julia toolkit.
__RELEASE__ -> Released date of this julia toolkit.
__AUTHORS__ -> Authors of this julia toolkit.
#
authors     -> Print the authors of NevanlinnaAC to screen.
```
=#

#
include("global.jl")
#
export Ftype, Ctype
#
export __LIBNAME__
export __VERSION__
export __RELEASE__
export __AUTHORS__
#
export authors


#=
### *Includes And Exports* : *types.jl*
=#

#=
*Summary* :

Define some dicts and structs, which are used to store the config
parameters or represent some essential data structures.

*Members* :

```text
OperatorType    -> Enumerated type for Fermionic/Bosonic correlators.
Fermi, Bose     -> Values of OperatorType
Options         -> Options of Nevanlinna AC algorithm
#
AbstractData    -> Abstract input data in imaginary axis.
RawData         -> Raw input data.
GenSchurData    -> Preprocessed input data that used in generalized schur algorithm.
NevData         -> Preprocessed input data that used in Nevanlinna algorithm.
```
=#

#
include("types.jl")
#
export OperatorType, Bose, Fermi
export Options
#
export RawData
export GenSchurData 
export NevData 


#=
### *Includes And Exports* : *math.jl*
=#

#=
*Summary* :

Provide conformal transforms and basis functions in the complex plane.

*Members* :

```text
eye    -> To construct identity matrix
#
lft    -> The linear_fractional_transform.
mt     -> The mobius transform.
imt    -> The inverse mobius transform
#
_mti   -> The mobius transform with center 1.0Im
_imti   -> The inverse mobius transform with center 1.0Im
```
=#


#
include("math.jl")


#=
### *Includes And Exports* : *base.jl*
=#

#=
*Summary* :

To provide basic workflow for the users of the NevanlinnaAC toolkit.

*Members* :

```text
#
read_to_RawData -> Read the input data to RawData struct.
read_to_NevData -> Read the input data to NevData struct.
#
toNevData       -> transform RawData to NevData.
toGenSchurData -> transform RawData/NevData to toGenSchurData.
make_mesh       -> Generate mesh for the calculated spectrum.
```
=#

#
include("base.jl")
#
export read_to_RawData
export read_to_NevData
#
export toNevData
export toGenSchurData
export make_mesh



#=
### *Includes And Exports* : *interpolates.jl*
=#

#=
*Summary* :

Interpolation algorithms.

*Members* :

```text
pick_matrix  -> To construct Pick matrix
isvalid      -> Check validity of Data
issolvable   -> Check Solvability of given data
#
_recursion      -> Recursion relation within the generalized Schur algorithm.
_inv_recursion  -> Inverse recursion relation
coefficient     -> Calculate the coefficient of the recursion relation in generalized_schur algorithm.
#
schur_parameter   -> Calculate Schur parameters 
generalized_schur -> Generalized Schur algorithm
nevanlinna        -> Nevanlinna Interpolation
#
spectral_function -> calculate spectral function
```
=#


#
include("interpolates.jl")
#
export isvalid
export issolvable
#
export generalized_schur
export nevanlinna
#
export spectral_function




#
include("optimization.jl")
#
export loss
export fft_derivative
#
end