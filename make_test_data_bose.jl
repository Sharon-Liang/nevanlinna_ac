using Random; Random.seed!()
using DelimitedFiles, Printf, Plots
using HCubature, LinearAlgebra
using NevanlinnaAC: OperatorType, Bose, Fermi

function Masubara_freq(n::Int64, β::Real, type::OperatorType)
    type == Bose ? N = 2n : N = 2n + 1
    return N*π/β
end


"""
two_pole_model_GF(z, a1, a2, e1, e2)

    Green's Function for two pole model
"""
function two_pole_model_GF(z, a1=0.1, a2=0.7, e1=0.7, e2=2.5)
    return a1/(z^2 - e1^2) + a2/(z^2 - e2^2)
end


"""
_fermi_distribution(e, β, μ = 0.5)
    
    Fermi distribution function
"""
function _fermi_distribution(e, β = 1.0, μ = -0.5)
    den = exp((e-μ)*β) + 1
    return 1.0/den
end

"""
_energy_density(px, py)

    ϵ(px, py) = -2t[cos(px) + cos(py)]
"""
_energy_density(px, py) =  -0.5*(cos(px) + cos(py))


"""
hubbard_model_rpa()
    RPA result of Masubara frequency GF for doppend Hubbard model, energy unit D = 4t =1
"""
function _gf_kernel(z, qx, qy, px, py, β = 1.0, μ = -0.5)
    ep = _energy_density(px, py)
    epq = _energy_density(px+qx, py+qy)

    num = _fermi_distribution(ep, β, μ) - _fermi_distribution(epq, β, μ)
    den = ep - epq + z
    if abs(z) == 0. && abs(ep - epq) <= 1.e-10
        res = -β*exp((ep-μ)*β)/ _fermi_distribution(ep, β, μ)^2
    else
        res = num/den
    end
    return res
end

function hubbard_model_rpa(z, qx, qy)
    f = p -> _gf_kernel(z, qx, qy, p[1], p[2])/(2π)^2
    return hcubature(f, [-π, -π], [π, π])[1]
end


#Generate Test Data
Nf = 100
β = 1
ωnlist = [Masubara_freq(n, β, Bose) for n=0:Nf]
Giωnlist = map(z->hubbard_model_rpa(1.0im * z, π, π), ωnlist)

gfile = "./data/giwn_hubbard_M_point_noise.txt"
open(gfile, "w") do file
    write(file, "        ωn                   Re.G(iωn)           Im.G(iωn)      \n")
    write(file, "--------------------  --------------------  -------------------- \n")
    #writedlm(file, [ωnlist real.(Giωnlist) imag.(Giωnlist)])

    writedlm(file, [ωnlist real.(Giωnlist).+ 1.e-8 .* rand(length(Giωnlist)) imag.(Giωnlist).+ 1.e-8 .* rand(length(Giωnlist))])
end


ωlist = [i for i in range(-2π, 2π, length=500)]
η = 0.05
Alist = map(z->-2*imag(hubbard_model_rpa(z + 1.0im * η, π, π)), ωlist)


afile = "./data/Aw_hubbard_M_point.txt"
open(afile, "w") do file
    write(file, "        ωn                      A(ω)      \n")
    write(file, "--------------------  --------------------\n")
    writedlm(file, [ωlist Alist])
end



f, _, _ = multi_gaussian(3, μrange=[1,3])






