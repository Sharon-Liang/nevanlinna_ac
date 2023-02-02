"""
    OperatorType Bose Fermi

Operator type of physical observables
"""
@enum OperatorType Bose Fermi


"""
    Options

Options of NevanlinnaAC.

### Members

* precsion -> Precision of float numbers
* otype    -> Type of correlators
* ngrid    -> Number of Masubara frequencies used
* nmesh    -> Number of mesh points
* wmax     -> Right boundary (maximum value)
* wmin     -> Left boundary (minimum)
* ifrev    -> If to reverse the exrtrapolation order
* η        -> infinitesimal imaginary part 
"""
@with_kw struct Options
    precision :: Int = 128
    otype :: OperatorType
    ngrid :: Int = 20
    nmesh :: Int = 500
    wmax :: Float64 = 2π
    wmin :: Float64 = -2π
    torev :: Bool = true
    η :: Float64 = 0.05
end


#=
### *Customized Structs* : *Input Data*
=#

"""
    AbstractData

An abstract type representing the input data in for various exrtrapolation algorithms. 
"""
abstract type AbstractData end


"""
    RawData

It represent the raw input data. The datatype of raw data is complex.

### Members

* grid -> Raw iωₙ grid for the input data.
* value -> Raw input data C(iωₙ) = -G(iωₙ).

See also: [`GenSchurData`](@ref), [`NevData`](@ref).
"""
mutable struct RawData{T} <: AbstractData
    grid :: Vector{T}
    value :: Vector{T}
end


"""
    GenSchurData

It represent the input data for generalized schur algorithm. The datatype of GenSchurData is complex.

### Members

* grid  -> Transformed raw input data that live in the upper half plane.
* value -> Transformed raw input data that live within the unit circle in the complex plane.

See also: [`RawData`](@ref), [`NevData`](@ref).
"""
mutable struct GenSchurData{T} <: AbstractData
    grid :: Vector{T}
    value :: Vector{T}
end


"""
    NevData

It represent the input data for Nevanlinna algorithm. The datatype of NevData is complex.

### Members

* grid  -> Transformed raw input data that live in the upper half plane.
* value -> Transformed raw input data that live within the unit circle on the complex plane.

See also: [`RawData`](@ref), [`GenSchurData`](@ref).
"""
mutable struct NevData{T} <: AbstractData
    grid :: Vector{T}
    value :: Vector{T}
end



