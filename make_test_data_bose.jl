using Random; Random.seed!()
using DelimitedFiles, Printf
using HCubature, LinearAlgebra
using NevanlinnaAC: OperatorType, Bose, Fermi

function Masubara_freq(n::Int64, β::Real, type::OperatorType)
    type == Bose ? N = 2n : N = 2n + 1
    return N*π/β
end


"""
two_pole_model(z, a1, a2, e1, e2)

    Green's Function for two pole model
"""
function two_pole_model(z, a1=0.1, a2=0.7, e1=0.7, e2=2.5)
    return a1/(z^2 - e1^2) + a2/(z^2 - e2^2)
end


"""
_fermi_distribution(e, β, μ = 0.5)
    
    Fermi distribution function
"""
function _fermi_distribution(e, β, μ = 0.)
    den = exp((e-μ)*β) + 1
    return 1.0/den
end


"""
_rpa_kernel(z, e1, e2, β, μ)

    _rpa_kernel = [nF(e1) - nF(e2)]/[z + e1 - e2] where `nF(e)` is the _fermi_distribution function
"""
function _rpa_kernel(z, e1, e2, β, μ)
    if abs(z) == 0. && abs(e1 - e2) <= 1.e-10
        ϵ = (e1-μ)*β/2
        res = -β * sech(ϵ)^2 /4
    else
        num = _fermi_distribution(e1, β, μ) - _fermi_distribution(e2, β, μ)
        den = e1 - e2 + z
        res = num/den
    end
    return res
end


"""
_energy_density(px, py)

    ϵ(px, py) = -2t[cos(px) + cos(py)]
"""
_energy_density(px, py) =  -0.5*(cos(px) + cos(py))



"""
hubbard_model(z, qx, qy, β=1.0, μ=-0.5)
    RPA result of Masubara frequency GF for doppend Hubbard model, energy unit D = 4t =1
"""
function _hubbard_gf_kernel(z, qx, qy, px, py, β, μ)
    ep = _energy_density(px, py)
    epq = _energy_density(px+qx, py+qy)
    res = _rpa_kernel(z, ep, epq, β, μ)
    return convert(typeof(z), res)
end

function hubbard_model(z, qx, qy, β=1.0, μ=-0.5)
    f = p -> _hubbard_gf_kernel(z, qx, qy, p[1], p[2], β, μ)/(2π)^2
    return hcubature(f, [-π, -π], [π, π])[1]
end


"""
bandgap_model(z, qx, qy, β=1.0, μ=1.7, eshift=3.)
    RPA result of Masubara frequency GF for bandgap_model, energy unit D = 4t =1
"""
function _band_gap_gf_kernel(z, qx, qy, px, py, β, μ, eshift)
    ep = _energy_density(px, py) 
    epq = _energy_density(px+qx, py+qy)
    res = 0.
    for e1 in [ep, ep + eshift], e2 in [epq, epq + eshift]
        res += _rpa_kernel(z, e1, e2, β, μ)
    end
    return res
end

function bandgap_model(z, qx, qy, β=1.0, μ=1.7, eshift=3.)
    f = p -> _band_gap_gf_kernel(z, qx, qy, p[1], p[2], β, μ, eshift)/(2π)^2
    return hcubature(f, [-π, -π], [π, π])[1]
end

