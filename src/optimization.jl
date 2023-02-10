Zygote.@nograd pick_matrix
Zygote.@nograd issolvable
Zygote.@nograd schur_parameter
Zygote.@nograd coefficient

Zygote.@nograd toNevData
Zygote.@nograd toGenSchurData
#Zygote.@nograd fftfreq
Zygote.@nograd make_mesh

"""
    gradient_function(loss, pars::AbstractArray)

Generate the gradient function of `loss`
"""
function gradient_function(loss, pars::AbstractArray)
    g! = function (g, pars)
        grads = Zygote.gradient(loss, pars)[1]
        copy!(g, grads)
    end
    return g!
end


function loss(params::AbstractArray, d::RawData, option::Options; λ::Real=1.e-4)
    @unpack wmax, otype = option

    wmesh, Aw = spectral_function(option, d, params)
    L = 2*wmax
    Δω = L / lastindex(wmesh)

    ∂²Aw = map(i->(Aw[i+2]+Aw[i]-2*Aw[i+1])/Δω^2, 1:lastindex(wmesh)-2) 
    smooth_condition = λ * norm(∂²Aw)^2 * Δω 
    
    if otype == Fermi
        sum_rule = abs(1.0 - sum(Aw) * Δω)^2
    else
        Ãw = @. Aw / wmesh
        pos = findall(x -> x==0.0, wmesh)
        for i in pos
            Ãw[i] = (Ãw[i-1] + Ãw[i+1])/2.0
        end

        ind = findall(x -> x==0.0, d.grid)[1]
        sum_rule = abs(-d.value[ind] - sum(Ãw) * Δω)^2
    end

    return smooth_condition + sum_rule
end
