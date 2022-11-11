using Random; Random.seed!()
using DelimitedFiles, Printf, Plots
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
function _hubbard_gf_kernel(z, qx, qy, px, py, β = 1.0, μ = -0.5)
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

function hubbard_model_rpa(z, qx, qy, β=1.0, μ=-0.5)
    f = p -> _hubbard_gf_kernel(z, qx, qy, p[1], p[2], β, μ)/(2π)^2
    return hcubature(f, [-π, -π], [π, π])[1]
end


"""
band_gap_model()
    RPA result of Masubara frequency GF for band_gap_model, energy unit D = 4t =1
"""
function _band_gap_gf_kernel(z, qx, qy, px, py, β = 1.0, μ = 1.7)
    ep_val = _energy_density(px, py) 
    epq_val = _energy_density(px+qx, py+qy)
    res = 0.
    for ep in [ep_val, ep_val + 3.], epq in [epq_val, epq_val + 3.]
        if abs(z) == 0. && abs(ep - epq) <= 1.e-10
            res += -β*exp((ep-μ)*β)/ _fermi_distribution(ep, β, μ)^2
        else
            num = _fermi_distribution(ep, β, μ) - _fermi_distribution(epq, β, μ)
            den = ep - epq + z
            res += num/den
        end
    end
    return res
end

function band_gap_model(z, qx, qy, β=1.0, μ=1.7)
    f = p -> _band_gap_gf_kernel(z, qx, qy, p[1], p[2], β, μ)/(2π)^2
    return hcubature(f, [-π, -π], [π, π])[1]
end

#Generate Test Data
β = 1

p1 = 0
p2 = π/2

ωlist = [i for i in range(-2π, 2π, length=500)]
η = 0.05
Alist = map(z->-2*imag(two_pole_model(z + 1.0im * η)), ωlist)

plot(ωlist, Alist)



afile = @sprintf "./data/Aw_two_pole_model_beta_%.2f.txt" β
open(afile, "w") do file
    write(file, "        ω                       A(ω)      \n")
    write(file, "--------------------  --------------------\n")
    writedlm(file, [ωlist Alist])
end


Nf = 100
ωnlist = [Masubara_freq(n, β, Bose) for n=0:Nf]
Giωnlist = map(z->two_pole_model(1.0im * z), ωnlist)

gfile = @sprintf "./data/giwn_two_pole_model_beta_%.2f.txt" β
open(gfile, "w") do file
    write(file, "        ωn                   Re.G(iωn)           Im.G(iωn)      \n")
    write(file, "--------------------  --------------------  -------------------- \n")
    writedlm(file, [ωnlist real.(Giωnlist) imag.(Giωnlist)])

    #writedlm(file, [ωnlist real.(Giωnlist).+ 1.e-8 .* rand(length(Giωnlist)) imag.(Giωnlist).+ 1.e-8 .* rand(length(Giωnlist))])
end









