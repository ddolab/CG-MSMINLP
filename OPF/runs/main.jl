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

@everywhere using OPF
using Random
using Distributions
using XLSX

include("instance.jl")
include("write_functions.jl")

model_type = ARGS[1]
idx = ARGS[2] # instance number

ITER_LIM = 1000
TIME_LIM = 10800.0
GAP = 1e-3
SP_GAP = 1e-4

seeds = Dict("1"=>1215,"2"=>1222,"3"=>1316,"4"=>1541,"5"=>1950,
    "6"=>1361,"7"=>1572,"8"=>1188,"9"=>1399,"10"=>1257,
    "11"=>2711,"12"=>2402,"13"=>2913,"14"=>2438,"15"=>2445,
    "16"=>2409,"17"=>2042,"18"=>2346,"19"=>4330,"20"=>4732,
    "21"=>2698,"22"=>2147,"23"=>2688,"24"=>2029,"25"=>2130,
    "26"=>3118,"27"=>3788,"28"=>3454,"29"=>3534,"30"=>3235,
    "31"=>3236,"32"=>3756,"33"=>3809,"34"=>3149,"35"=>3040,
    "36"=>3234,"37"=>3586,"38"=>3390,"39"=>3930,"40"=>3191)

function instance_specifications(instance_number::Int)
    num_gensets, num_periods, num_hours, br_structure = 0, 0, 24, []

    if instance_number in 1:20
        num_gensets = 1
        if instance_number in 1:5
            num_periods = 2
            br_structure = [1,3,3]
        elseif instance_number in 6:10
            num_periods = 3
            br_structure = [1,3,3,4]
        elseif instance_number in 11:15
            num_periods = 4
            br_structure = [1,3,3,3,3]
        elseif instance_number in 16:20
            num_periods = 5
            br_structure = [1,4,4,2,2,2]
        end
    elseif instance_number in 21:40
        num_gensets = 2
        if instance_number in 21:25
            num_periods = 2
            br_structure = [1,3,3]
        elseif instance_number in 26:30
            num_periods = 3
            br_structure = [1,3,3,4]
        elseif instance_number in 31:35
            num_periods = 4
            br_structure = [1,3,3,3,3]
        elseif instance_number in 36:40
            num_periods = 5
            br_structure = [1,4,4,2,2,2]
        end
    else
        error("Invalid instance number")
    end 

    return num_gensets, num_periods, num_hours, br_structure
end

nx = generate_network()
nx_data = network_data(nx.buses, nx.lines, "network_data.xlsx")

num_gensets, num_periods, num_hours, br_structure = instance_specifications(parse(Int, idx))
tree = Tree(br_structure)

data = generate_data(nx, nx_data, num_gensets, num_periods, num_hours, br_structure; seed = seeds[idx])

sol = 0
if model_type == "1"      #Fullspace model
    sol = solve_stoch_model_FS(data; tree = tree, timelimit = TIME_LIM, gap = GAP)
elseif model_type == "2"  #CG
    sol = solve_CG(data; tree = tree, iterlimit = ITER_LIM, timelimit = TIME_LIM, relative_gap = GAP, sp_gap = SP_GAP, do_column_sharing = false)
elseif model_type == "3"  #CG with column sharing
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