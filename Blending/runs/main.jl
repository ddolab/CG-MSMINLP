using Distributed

if ARGS[1] == "1" # Fullspace model
    using Pkg
    if isfile("Project.toml") && isfile("Manifest.toml")
        Pkg.activate(".")
    end
else # CG/CGCS
    if nworkers() > 1
        rmprocs(workers())
    end
    addprocs(64)

    @everywhere begin
        using Pkg
        if isfile("Project.toml") && isfile("Manifest.toml")
            Pkg.activate(".")
        end
    end
end

@everywhere using Blending
using Random
using Distributions
using XLSX

include("instance.jl")
include("write_functions.jl")

model_type = ARGS[1]
idx = ARGS[2] # instance number

ITER_LIM = 1000
TIME_LIM = 7200.0
GAP = 1e-4
SP_GAP = 1e-4

seeds = Dict("1"=>12345, "2"=>12346, "3"=>12347, "4"=>12348, "5"=>12349, 
    "6"=>123410,"7"=>123411, "8"=>123412, "9"=>123413, "10"=>123414, 
    "11"=>123415, "12"=>123416, "13"=>123417, "14"=>123418, "15"=>123419,
    "16"=>923435, "17"=>923436, "18"=>923437, "19"=>923438, "20"=>923439,
    "21"=>123420, "22"=>123421, "23"=>123422, "24"=>123423, "25"=>123424,
    "26"=>123425, "27"=>123426, "28"=>123427, "29"=>123428, "30"=>123429,
    "31"=>1234301, "32"=>123431, "33"=>123432, "34"=>123433, "35"=>123434,
    "36"=>923440, "37"=>923441, "38"=>923442, "39"=>923443, "40"=>923446)

# This function is used to generate the branching structure for the tree
function branching_structure(num_periods::Int)
    br_structure = [1]
    for i in 1:num_periods
        push!(br_structure, 2)
    end
    return br_structure
end

function instance_data(instance_number::Int)
    num_inputs, num_blends, num_periods = 0, 0, 0

    if instance_number in 1:20
        num_inputs = 5
        num_blends = 3
        if instance_number in 1:5
            num_periods = 3
        elseif instance_number in 6:10
            num_periods = 5
        elseif instance_number in 11:15
            num_periods = 7
        elseif instance_number in 16:20
            num_periods = 9
        end
    elseif instance_number in 21:40
        num_inputs = 12
        num_blends = 10
        if instance_number in 21:25
            num_periods = 3
        elseif instance_number in 26:30
            num_periods = 5
        elseif instance_number in 31:35
            num_periods = 7
        elseif instance_number in 36:40
            num_periods = 9
        end
    else
        error("Invalid instance number")
    end 

    return [num_inputs, num_blends, num_periods, branching_structure(num_periods)]
end

num_inputs, num_blends, num_periods, br_structure = instance_data(parse(Int, idx)) 
@assert length(br_structure) == num_periods+1

tree = Tree(br_structure)

data = generate_data(num_inputs, num_blends, num_periods, br_structure; seed = seeds[idx])


sol = 0
if model_type == "1" # Fullspace model
    sol = solve_stoch_model_FS(data; tree = tree, timelimit = TIME_LIM, gap = GAP)
elseif model_type == "2" # CG model
    sol = solve_CG(data; tree = tree, iterlimit = ITER_LIM, timelimit = TIME_LIM, relative_gap = GAP, sp_gap = SP_GAP, do_column_sharing = false)
elseif model_type == "3" # CGCS model
    sol = solve_CG(data; tree = tree, iterlimit = ITER_LIM, timelimit = TIME_LIM, relative_gap = GAP, sp_gap = SP_GAP, do_column_sharing = true)
else
    error("Invalid model type")
end

if model_type == "1"
    create_spreadsheet_FS(sol, model_type, idx)
else
    rmprocs(workers())
    create_spreadsheet_CG(sol, model_type, idx)
    write_CG_gap_to_file(sol, model_type, idx)
end

