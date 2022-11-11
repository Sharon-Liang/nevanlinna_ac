using DelimitedFiles, Printf, LinearAlgebra, Random; Random.seed!()
using NevanlinnaAC, Optim, Zygote
using NevanlinnaAC: loss_fermi, gradient_function
using Dates, HDF5, Printf, ArgParse
include("optim_outputs.jl")

function read_GF_data(path::String)
    data = readdlm(path)
    x = map(x->convert(eltype(data), x), data[2:end,1]) 
    y = data[:,2] .+ 1.0im .* data[:,3]
    return x, y
end

settings = ArgParseSettings(prog="loss fermi"
)
@add_arg_table! settings begin
    "--prec"
        arg_type = Int64
        default = 128
        help = "Float number precision"
    "--ωmax"
        arg_type = Float64
        default = 4π
        help = "maximum ω value"    
    "--Nω"
        arg_type = Int64
        default = 1000
        help = "number of ω"     
    "--eta"
        arg_type = Float64
        default = 0.05
        help = "ω + iη"     
    "--lambda"
        arg_type = Float64
        default = 1.e-4
        help = "prefactor on smooth_condition"  
    "--Nf"
        arg_type = Int64
        default = 20
        help = "number of Masubara frequencies"    
    "--H"
        arg_type = Int64
        default = 20
        help = "order of hardy basis"  
    "--model"
        arg_type = String
        default = "gaussian_1"
        help = "data file name = giwn_model"  
    "--tag"
        arg_type = String
        default = Dates.format(now(), "yyyy-mm-dd")
        help = "date tag"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const prec = parsed_args[:prec]
setprecision(BigFloat, prec)
const float_type = BigFloat

const ωmax = parsed_args[:ωmax]
const Nω = parsed_args[:Nω]
const η = parsed_args[:eta]
const λ = parsed_args[:lambda]

const Nf = parsed_args[:Nf]
const H = parsed_args[:H]
const tag = parsed_args[:tag]

const model = parsed_args[:model]


DataFilePath = "/home/sliang/JuliaCode/NevanlinnaAC/data/gaussian/giwn_$(model).txt"
fermiX, fermiY = read_GF_data(DataFilePath)

ResultFolder = @sprintf "/data/sliang/NevanlinnaAC/%s" model
isdir(ResultFolder) || mkdir(ResultFolder)
Newloc = "$(ResultFolder)/giwn_$(model).txt"
isfile(Newloc) || cp(DataFilePath, "$(ResultFolder)/giwn_$(model).txt")

file_name = @sprintf "Nf_%03i_H_%02i_prec_%i_lambda_%.2e_eta_%.3f_Nw_%i" Nf H prec λ η Nω

params_file = @sprintf "%s/optim_params_%s.hdf5" ResultFolder file_name

loss = params -> loss_fermi(params, fermiX[1:Nf], fermiY[1:Nf]; ωmax, Nω, η, λ, show_warning=true)

p0 = zeros(2H)
maxiter = 2000
save_every = 5
key = string(0)
h5open(params_file, "w") do file 
    g = create_group(file, key)
    g["params"] = p0
    g["grads"] = Zygote.gradient(loss, p0)[1]
end

for itr = 1: div(maxiter, save_every)
    res = Optim.optimize(loss, gradient_function(loss, p0), p0, LBFGS(), Optim.Options(show_trace=true, store_trace=true, iterations = save_every))

    if Optim.converged(res) == true
        global p0 = res.minimizer

        h5open(params_file, "r+") do file 
            key = (itr - 1) * save_every + Optim.iterations(res) |> Int |> string
            g = create_group(file, key)
            g["params"] = p0
            g["grads"] = Zygote.gradient(loss, p0)[1]
        end
        break
    else
        global p0 = res.minimizer
        h5open(params_file, "r+") do file 
            key = itr * save_every |> Int |> string
            g = create_group(file, key)
            g["params"] = p0
            g["grads"] = Zygote.gradient(loss, p0)[1]
        end
    end
end

w0, A0 = spectral_function(Fermi, fermiX[1:Nf], fermiY[1:Nf]; ωmax, Nω, η, float_type);
w1, A1 = spectral_function(Fermi, fermiX[1:Nf], fermiY[1:Nf], p0; ωmax, Nω, η, float_type);

spectral_function_file = @sprintf "%s/Aw_%s.txt" ResultFolder file_name
open(spectral_function_file, "w") do file
    write(file, "        ω                   A(ω)-inite-0      A(ω)-hardy-optim  \n")
    write(file, "--------------------  --------------------  --------------------\n")
    writedlm(file, [w0 A0 A1])
end