using Plots
include("./make_test_data_bose.jl")

#G for Γ
@enum Point M X G
function momentum(p::Point)
    if p == M
        px = π
        py = π
    elseif p == X
        px = π
        py = π/2
    else
        px = 0.
        py = 0
    end
    return px, py
end

η = 0.05;
ωlist = [i for i in range(-2π, 2π, length=500)];
Nf = 100;

Plist = [X, M, G]
betalist = [1.0, 10.0, 100.0]

modellist = [:hubbard, :bandgap]
for model in modellist
    if model ==:hubbard
        func = hubbard_model
    elseif model == :bandgap 
        func = bandgap_model
    end
    
    println("$(model)")
    for β in betalist, P in Plist
        println("beta=$(β), point=$(P)")

        px, py = momentum(P)
        Alist = map(z->-2*imag(func(z + 1.0im * η, px, py, β)), ωlist)

        ωnlist = [Masubara_freq(n, β, Bose) for n=0:Nf]
        Giωnlist = map(z -> func(1.0im * z, px, py, β), ωnlist)

        afile = @sprintf "./data/Aw_%s_%s_point_beta_%.2f.txt" model P β
        open(afile, "w") do file
            write(file, "        ω                       A(ω)      \n")
            write(file, "--------------------  --------------------\n")
        writedlm(file, [ωlist Alist])
        end
        gfile = @sprintf "./data/giwn_%s_%s_point_beta_%.2f.txt" model P β
        open(gfile, "w") do file
            write(file, "        ωn                   Re.G(iωn)           Im.G(iωn)      \n")
            write(file, "--------------------  --------------------  -------------------- \n")
            writedlm(file, [ωnlist real.(Giωnlist) imag.(Giωnlist)])
        end
    end
end