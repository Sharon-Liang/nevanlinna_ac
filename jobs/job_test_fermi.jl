using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

tag = Dates.format(now(), "yyyy-mm-dd")

env = "/home/sliang/JuliaCode/NevanlinnaAC"
prog = "/home/sliang/JuliaCode/NevanlinnaAC/jobs/test_fermi.jl"

ωmax = 4π
Nωlist = [500, 1000]
ηlist = [0.05]
λlist = [10^i for i in range(-8, -1, step=1)]
Nflist = [10, 20, 30, 40]
Hlist = [i for i in range(5, 30, step=5)]

modellist = ["gaussian_1", "gaussian_3", "delta_eta_0.05"]
preclist = [128, 256]


for model in modellist
    ResultFolder = @sprintf "/data/sliang/NevanlinnaAC/%s" model
    isdir(ResultFolder) || mkdir(ResultFolder)
    logFolder = "$(ResultFolder)/log"
    isdir(logFolder) || mkdir(logFolder)

    for prec in preclist, Nω in Nωlist, eta in ηlist, lambda in λlist, Nf in Nflist, H in Hlist
        args = Dict("model" => model,
                "prec" => prec,
                "ωmax" => ωmax,
                "Nω" => Nω,
                "eta" => eta,
                "lambda" => lambda, 
                "Nf" => Nf,
                "H" => H,
                "tag" => tag
            )
        file_name = @sprintf "Nf_%03i_H_%02i_prec_%i_lambda_%.2e_eta_%.3f_Nw_%i" Nf H prec lambda eta Nω

        logdir =  logFolder * "/" * file_name * "_" * tag
        submitJob(env, prog, args, logdir, Run = true, partitian = v100)
    end
end