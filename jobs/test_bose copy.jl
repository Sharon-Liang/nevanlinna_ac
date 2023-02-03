using LinearAlgebra; BLAS.set_num_threads(Threads.nthreads())
using TimerOutputs, Dates, ArgParse
using DelimitedFiles, HDF5, Printf
using Random; Random.seed!()

using NevanlinnaAC, Optim, Zygote
using NevanlinnaAC: loss_bose, gradient_function
include("optim_outputs.jl")


const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

function read_GF_data(path::String)
    data = readdlm(path, skipstart = 2)
    x = map(x->convert(ComplexF64, x), data[2:end,1]) 
    y = data[2:end,2] .+ 1.0im .* data[2:end,3]
    g0 = data[1,2]
    return g0, x, y
end

settings = ArgParseSettings(prog="loss bose"
)
@add_arg_table! settings begin
    "--prec"
        arg_type = Int64
        default = 128
        help = "Float number precision"
    "--use_g0"
        arg_type = Bool
        default = false
        help = "whether to G(iωn=0) = g0 or not"    
    "--isodd"
        arg_type = Bool
        default = true
        help = "whether to use the odd_condition"
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
        default = "hubbard_M_point"
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

const isodd = parsed_args[:isodd]
const use_g0 = parsed_args[:use_g0]

const ωmax = parsed_args[:ωmax]
const Nω = parsed_args[:Nω]
const η = parsed_args[:eta]
const λ = parsed_args[:lambda]

const Nf = parsed_args[:Nf]
const H = parsed_args[:H]
const tag = parsed_args[:tag]

const model = parsed_args[:model]


DataFilePath = "/home/sliang/JuliaCode/NevanlinnaAC/data/giwn_$(model).txt"
g0, boseX, boseY = read_GF_data(DataFilePath)

ResultFolder = @sprintf "/data/sliang/NevanlinnaAC/%s" model
isdir(ResultFolder) || mkdir(ResultFolder)
giwn_Newloc = "$(ResultFolder)/giwn_$(model).txt"
isfile(giwn_Newloc) || cp(DataFilePath, "$(ResultFolder)/giwn_$(model).txt")

Aw_DataFilePath =  "/home/sliang/JuliaCode/NevanlinnaAC/data/Aw_$(model).txt"
Aw_Newloc = "$(ResultFolder)/Aw_$(model).txt"
isfile(Aw_Newloc) || cp(Aw_DataFilePath, "$(ResultFolder)/Aw_$(model).txt")

file_name = @sprintf "g0_%s_odd_%s_Nf_%03i_H_%02i_prec_%i_lambda_%.2e_eta_%.2f_Nw_%i" use_g0 isodd Nf H prec λ η Nω

params_file = @sprintf "%s/optim_params_%s.hdf5" ResultFolder file_name


loss = params -> loss_bose(params, boseX[1:Nf], boseY[1:Nf], g0; ωmax, Nω, η, λ, isodd, use_g0, show_warning=true)

p0 = zeros(2H)
maxiter = 1200
save_every = 5
key = string(0)
h5open(params_file, "w") do file 
    g = create_group(file, key)
    g["params"] = p0
    g["grads"] = Zygote.gradient(loss, p0)[1]
end

for itr = 1: div(maxiter, save_every)
    @timeit to "optimize" begin
        res = Optim.optimize(loss, gradient_function(loss, p0), p0, LBFGS(), Optim.Options(show_trace=true, store_trace=true, iterations = save_every))

        if Optim.converged(res) == true
            global p0 = res.minimizer

            h5open(params_file, "r+") do file 
                key = (itr - 1) * save_every + Optim.iterations(res) |> Int |> string
                g = create_group(file, key)
                g["params"] = p0
                g["grads"] = Zygote.gradient(loss, p0)[1]
            end

            println(res)
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
end



w0, A0 = spectral_function(Bose, boseX[1:Nf], boseY[1:Nf], g0; ωmax, Nω, η, float_type, use_g0);
w1, A1 = spectral_function(Bose, boseX[1:Nf], boseY[1:Nf], g0, p0; ωmax, Nω, η, float_type, use_g0);

spectral_function_file = @sprintf "%s/Aw_%s.txt" ResultFolder file_name
open(spectral_function_file, "w") do file
    write(file, "        ω                   A(ω)-inite-0      A(ω)-hardy-optim  \n")
    write(file, "--------------------  --------------------  --------------------\n")
    writedlm(file, [w0 A0 A1])
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")