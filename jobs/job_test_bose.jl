using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

tag = Dates.format(now(), "yyyy-mm-dd")

env = "/home/sliang/JuliaCode/NevanlinnaAC"
prog = "/home/sliang/JuliaCode/NevanlinnaAC/jobs/test_bose.jl"

ωmax = 4π
Nωlist = [500, 1000]
ηlist = [0.05]
λlist = [10^i for i in range(-8, -1, step=1)]
Nflist = [10, 20, 30, 40]
Hlist = [i for i in range(5, 30, step=5)]

modellist = ["hubbard_M_point_beta_1.00", "bandgap_M_point_beta_10.00", "bandgap_X_point_beta_10.00"]
preclist = [128, 256, 512]

use_g0_list = [false, true]
isodd_list = [false, true]

#====test params====#
#Nωlist = [500]
#ηlist = [0.05]
#λlist = [1.e-4]
#Nflist = [10]
#Hlist = [5]

#modellist = ["hubbard_M_point_beta_1.00"]
#preclist = [128]

#use_g0_list = [false]
#isodd_list = [false]

for model in modellist
    ResultFolder = @sprintf "/data/sliang/NevanlinnaAC/%s" model
    isdir(ResultFolder) || mkdir(ResultFolder)
    logFolder = "$(ResultFolder)/log"
    isdir(logFolder) || mkdir(logFolder)

    for prec in preclist, Nω in Nωlist, eta in ηlist, lambda in λlist, Nf in Nflist, H in Hlist, use_g0 in use_g0_list, isodd in isodd_list
        if use_g0 == true || isodd == true
            args = Dict("model" => model,
                    "prec" => prec,
                    "isodd" => isodd,
                    "use_g0" => use_g0,
                    "ωmax" => ωmax,
                    "Nω" => Nω,
                    "eta" => eta,
                    "lambda" => lambda, 
                    "Nf" => Nf,
                    "H" => H,
                    "tag" => tag
                )
            file_name = @sprintf "%s_g0_%s_odd_%s_Nf_%03i_H_%02i_prec_%i_lambda_%.2e_eta_%.2f_Nw_%i" model use_g0 isodd Nf H prec lambda eta Nω

            logdir =  logFolder * "/" * file_name * "_" * tag
            submitJob(env, prog, args, logdir, Run = true, partitian = v100)
        end
    end
end