function generate_network()
    buses = 1:15 # bus indices starts from 1 instead of 0
    lines = 2:15
    children = Dict(1 => [2,13],
                    2 => [3],
                    3 => [4],
                    4 => [5,9],
                    5 => [6],
                    6 => [7],
                    7 => [],
                    8 => [],
                    9 => [8,10],
                    10 => [11],
                    11 => [12],
                    12 => [],
                    13 => [14],
                    14 => [15],
                    15 => [])
    parent = Dict(1 => [],
                  2 => [1],
                  3 => [2],
                  4 => [3],
                  5 => [4],
                  6 => [5],
                  7 => [6],
                  8 => [9],
                  9 => [4],
                  10 => [9],
                  11 => [10],
                  12 => [11],
                  13 => [1],
                  14 => [13],
                  15 => [14])
    return TestNetwork(buses, lines, parent, children)
end

function network_data(buses, lines, filename)
    xf = XLSX.readxlsx(filename)
    sh = xf["Data"]
    r, x, A, pd_base, qd_base, b, g, Cp_base, C_VoLL, pg_base, qg_base, v_lb, v_ub = [[] for _ = 1:13]
    cell_ref = Dict("r"=>"B",
        "x"=> "C",
        "A"=> "D",
        "pd_base"=> "E",
        "qd_base"=> "F",
        "b"=> "G",
        "g"=> "H",
        "Cp_base"=> "I",
        "C_VoLL"=> "J",
        "pg_base"=> "K",
        "qg_base"=> "L",
        "v_lb"=> "M",
        "v_ub"=> "N" )
    for (key, value) in cell_ref
        for bus in buses
            if key == "r"
                push!(r, sh[string(value, bus+1)])
            elseif key == "x"
                push!(x, sh[string(value, bus+1)])
            elseif key == "A"
                push!(A, sh[string(value, bus+1)])
            elseif key == "pd_base"
                push!(pd_base, sh[string(value, bus+1)])
            elseif key == "qd_base"
                push!(qd_base, sh[string(value, bus+1)])
            elseif key == "b"
                push!(b, sh[string(value, bus+1)])
            elseif key == "g"
                push!(g, sh[string(value, bus+1)])
            elseif key == "Cp_base"
                push!(Cp_base, sh[string(value, bus+1)])
            elseif key == "C_VoLL"
                push!(C_VoLL, sh[string(value, bus+1)])
            elseif key == "pg_base"
                push!(pg_base, sh[string(value, bus+1)])
            elseif key == "qg_base"
                push!(qg_base, sh[string(value, bus+1)])
            elseif key == "v_lb"
                push!(v_lb, sh[string(value, bus+1)])
            elseif key == "v_ub"
                push!(v_ub, sh[string(value, bus+1)])
            else
                error("$(key) does not exist!")
            end
        end
    end

    return TestSystemData(buses,lines,r,x,A,pd_base,qd_base,b,g,Cp_base,C_VoLL, pg_base,qg_base,v_lb,v_ub)
end


function generate_data(network::TestNetwork, network_data::TestSystemData, num_gensets::Int64, num_periods::Int64, num_hours::Int64, brstructure::Vector{Int64}; seed::Int64)

    @assert network.buses == network_data.buses
    @assert network.lines == network_data.lines

    rng = MersenneTwister(seed)

    total_buses = length(network_data.buses)
    
    capacity_of_genset = sort(round.(rand(rng, Uniform(0.1, 0.8), num_gensets); digits=5))

    transport_cost = round.(rand(rng, Uniform(25, 35), (num_gensets, num_periods)); digits=5)
    sort!(transport_cost, dims=1)

    genset_fuel_cost = round.(rand(rng, Uniform(100, 120), num_periods); digits=5)

    #--------------------------------------------------------------#
    time_of_day_factor = rand(rng, Uniform(0.5, 1), (num_periods, num_hours))
    # normalizing time of day factor s.t. atleast one 'h' has peak demand
    for t in 1:num_periods
        max_tod_factor = maximum(time_of_day_factor[t,:])
        for h in 1:num_hours
            time_of_day_factor[t,h] = time_of_day_factor[t,h]/max_tod_factor
        end
    end
    #--------------------------------------------------------------#

    #--------------------------------------------------------------#
    power_generation_cost = rand(rng, Uniform(0.8, 1.2), (total_buses, num_periods, num_hours))
    for i in 1:total_buses
        for t in 1:num_periods
            for h in 1:num_hours
                power_generation_cost[i,t,h] *= network_data.Cp_base[i]
            end

            #Making cost of power generation proportional to the time of day factor
            idx = sortperm(time_of_day_factor[t,:])
            new_power_generation_cost = similar(power_generation_cost[i,t,:])
            sorted_power_generation_cost = sort(power_generation_cost[i,t,:])
            for h in 1:num_hours
                new_power_generation_cost[idx[h]] = sorted_power_generation_cost[h]
            end
            power_generation_cost[i,t,:] = new_power_generation_cost
        end
    end
    power_generation_cost = round.(power_generation_cost; digits=5)
    #--------------------------------------------------------------#

    #--------------------------------------------------------------#
    active_pg_ub = rand(rng, Uniform(0.8,1.2), total_buses)
    reactive_pg_ub = rand(rng, Uniform(0.8,1.2), total_buses)
    for i in 1:total_buses
        active_pg_ub[i] *= network_data.pg_base[i]
        reactive_pg_ub[i] *= network_data.qg_base[i]
    end
    active_pg_ub = round.(active_pg_ub; digits=5)
    reactive_pg_ub = round.(reactive_pg_ub; digits=5)
    #--------------------------------------------------------------#

    alpha = ones(num_periods).*(168/num_hours)

    #------------------stochastic power demand data------------------#
    tree = Tree(brstructure)

    load_variation_factor = Vector{Float64}[]

    for n in 1:maximum(OPF.nodes(tree))
        if n == 1
            push!(load_variation_factor, Float64[])
        else
            push!(load_variation_factor, rand(rng, Uniform(0.5,1.5), total_buses))
        end
    end

    active_power_demand = Matrix{Float64}[]
    reactive_power_demand = Matrix{Float64}[]
    for n in 1:maximum(OPF.nodes(tree))
        push!(active_power_demand, Matrix{Float64}(undef, total_buses, num_hours))
        push!(reactive_power_demand, Matrix{Float64}(undef, total_buses, num_hours))
        if n == 1
            continue
        end

        for i in 1:total_buses
            for h in 1:num_hours
                active_power_demand[n][i,h] = round(network_data.pd_base[i]*load_variation_factor[n][i]*time_of_day_factor[OPF.stage(tree,n)[1],h]; digits=5)
                reactive_power_demand[n][i,h] = round(network_data.qd_base[i]*load_variation_factor[n][i]*time_of_day_factor[OPF.stage(tree,n)[1],h]; digits=5)
            end
        end
    end
    #--------------------------------------------------------------#

    # updating proabalities generated by the original Tree struct
    for n in OPF.nodes(tree)
        tree.probability[n] = 1/(brstructure[OPF.stage(tree, n)[1]+1])
    end

    node_total_prob = Float64[]
    for n in OPF.nodes(tree)
        push!(node_total_prob, round(OPF.total_node_probability(tree, n);digits=7))
    end

    return ModelDataNodalStochastic(network_data.buses,
        network_data.lines, 
        network.parent,
        network.children,
        num_gensets,
        num_periods,
        num_hours,
        alpha,
        transport_cost,
        power_generation_cost,
        network_data.C_VoLL,
        active_power_demand,
        reactive_power_demand,
        genset_fuel_cost,
        capacity_of_genset,
        network_data.r,
        network_data.g,
        network_data.x,
        network_data.A,
        network_data.b,
        active_pg_ub,
        reactive_pg_ub,
        network_data.v_lb,
        network_data.v_ub,
        node_total_prob)
end

