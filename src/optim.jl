function H2basis(z::Number, k::Int64)
    return 1/(√π*(z+1.0im)) * h(z)^k
end

function hardy_expand(z::Number, H::Int64, param::Matrix)
    res = 0.
    for k = 1:H
        Bk = H2basis(z, k)
        res += param[k,1] * Bk + param[k,2] * Bk'
    end
    return res
end

function loss(Aω::Vector)
    
end


# 2nd order derivative