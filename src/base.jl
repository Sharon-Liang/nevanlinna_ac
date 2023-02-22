#=
### *Transform between data types*
=#


"""
    toNevData(d::RawData, option::Options)

Convert RawData to NevData. 

* For fermionic Green's function: G(iωₙ) -> -G(iωₙ)
* For Bosonic Green's function: G(iωₙ) -> -iωₙ*G(iωₙ)
"""
function toNevData(d::RawData, option::Options)
    @unpack otype, η = option
    
    grid = copy(d.grid)

    if otype == Fermi
        value = -d.value
    else
        #When n=0, replace iωₙ by iη
        ind = findall(x -> x==0.0, grid)[1]
        grid[ind] = one(Ctype)im * Ctype(η) 
        value = - grid .* d.value
    end

    return NevData(grid, value)
end


"""
    toGenSchurData(d::RawData, option::Options)

Convert RawData to GenSchurData. 

* For fermionic Green's function: G(iωₙ) -> mti[-G(iωₙ)]
* For Bosonic Green's function: G(iωₙ) -> mti[-iωₙ*G(iωₙ)]
"""
function toGenSchurData(d::RawData, option::Options)
    nd = toNevData(d, option)
    value = map(_mti, nd.value)
    return GenSchurData(nd.grid, value)
end


"""
    toGenSchurData(d::NevData)

Convert NevData to GenSchurData. ``d.value -> mti(d.value)``

*
"""
toGenSchurData(d::NevData) = GenSchurData(d.grid, map(_mti, d.value))


#=
### *Read data*
=#

"""
    read_to_RawData(path::String, option::Options; kwargs...)

Read data in imaginary axis and return a `RawData` struct. A raw data file should be consist of 3 columns, the first column if Masubara frequency ωₙ, the second and the third column is the real and imaginary part of the Green's fuction G(iωₙ) respectively. The Masubara frequency should be started from n=1 for fermions and n=0 for bosons.

See also: [`RawData`](@ref).
"""

function read_to_RawData(path::String, option::Options; kwargs...)
    @unpack ngrid, otype, torev, η = option

    d = readdlm(path; kwargs...)
    grid = d[1:ngrid, 1] * one(eltype(d))im
    value = d[1:ngrid, 2] + d[1:ngrid, 3] * one(eltype(d))im
    
    grid = convert.(Ctype, grid)
    value = convert.(Ctype, value)

    if torev
        grid = reverse(grid)
        value = reverse(value)
    end
    
    return RawData(grid, value)
end


"""
    read_to_NevData(path::String, option::Options; kwargs...)

Read data in imaginary axis and return a `NevData` struct. A raw data file should be consist of 3 columns, the first column if Masubara frequency ωₙ, the second and the third column is the real and imaginary part of the Green's fuction G(iωₙ) respectively. The Masubara frequency should be started from n=1 for fermions and n=0 for bosons.

See also: [`NevData`](@ref).
"""
function read_to_NevData(path::String, option::Options; kwargs...)
    d = read_to_RawData(path, option; kwargs...)
    
    return toNevData(d, option)
end


"""
    make_mesh(option::Options)

Make mesh of real frequencies, return ω+iη. For Bosonic spectral, ``nmesh`` is rest to be the nearst odd number.
"""
function make_mesh(option::Options)
    @unpack nmesh, wmax, η = option
    wmax = Ctype(wmax)
    mesh = [i for i in range(-wmax, wmax, length = nmesh)]

    @. mesh += Ctype(1.0im*η)
    return mesh
end