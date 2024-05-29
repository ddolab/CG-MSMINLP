# For generating initial columns for the column generation procedure
function deterministic_model(data::ModelDataNodalStochastic;
    tree,
    timelimit,
    gap,
    desired_num_sols)
    
    pd, qd = expected_power_demand(data; scenario_tree = tree)

    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "TimeLimit", timelimit)
    set_optimizer_attribute(model, "MIPGap", gap)
    set_optimizer_attribute(model, "PoolSolutions", desired_num_sols)
    set_optimizer_attribute(model, "PoolSearchMode", 1)

    # define variables
    @variables(
        model,
        begin
            # binary variables
            y[m in 1:data.M, i in data.B, t in 1:data.T], Bin
            z[m in 1:data.M, t in 2:data.T], Bin
            
            # (semi-)bounded continuous variables
            0 <= pg[i in data.B, t in 1:data.T, h in 1:data.H] <= data.pg_ub[i]
            0 <= qg[i in data.B, t in 1:data.T, h in 1:data.H] <= data.qg_ub[i]
            0 <= delta[i in data.B, t in 1:data.T, h in 1:data.H] <= 1
            data.v_lb[i] <= v[i in data.B, t in 1:data.T, h in 1:data.H] <= data.v_ub[i]
            l[i in data.L, t in 1:data.T, h in 1:data.H] >= 0
            p_a[m in 1:data.M, i in data.B, t in 1:data.T, h in 1:data.H] >= 0

            # unrestricted
            q_a[m in 1:data.M, i in data.B, t in 1:data.T, h in 1:data.H]
            P[i in data.L, t in 1:data.T, h in 1:data.H]
            Q[i in data.L, t in 1:data.T, h in 1:data.H]
        end
    )

    @constraint(
        model,
        [m in 1:data.M, t in 1:data.T], sum(y[m, i, t] for i in data.B)  == 1
    )

    @constraints(
        model,
        begin
            [m in 1:data.M, i in data.B, t in 2:data.T], y[m,i,t]-y[m,i,t-1] <= z[m,t]
            [m in 1:data.M, i in data.B, t in 2:data.T], z[m,t] <= 2 - y[m,i,t] - y[m,i,t-1] 
        end
    )

    # Power generation constraint
    @constraint(
        model,
        [m in 1:data.M, i in data.B, t in 1:data.T, h in 1:data.H],
        p_a[m,i,t,h] <= data.pa[m]*y[m,i,t]
    )
    @constraints(model, begin
        [m in 1:data.M, i in data.B, t in 1:data.T, h in 1:data.H], q_a[m,i,t,h] <= p_a[m,i,t,h] 
        [m in 1:data.M, i in data.B, t in 1:data.T, h in 1:data.H], q_a[m,i,t,h] >= -p_a[m,i,t,h]
    end)

    # Power flow equations
    @constraints(
        model,
        begin
            [i in data.B, t in 1:data.T, h in 1:data.H],
                pg[i,t,h] - pd[t][i,h]*delta[i,t,h] + sum(p_a[m,i,t,h] for m in 1:data.M) ==
                sum(P[j,t,h] for j in data.children[i]) - sum((P[i,t,h]-data.r[i]*l[i,t,h]) for _ in data.parent[i]) + data.g[i]*v[i,t,h] 
            
            [i in data.B, t in 1:data.T, h in 1:data.H],
            qg[i,t,h] - qd[t][i,h]*delta[i,t,h] + sum(q_a[m,i,t,h] for m in 1:data.M) ==
            sum(Q[j,t,h] for j in data.children[i]) - sum((Q[i,t,h]-data.x[i]*l[i,t,h]) for _ in data.parent[i]) + data.b[i]*v[i,t,h] 
        end
    )

    @constraints(
        model, 
        begin
            [i in data.L, t in 1:data.T, h in 1:data.H],
            v[i,t,h] == v[data.parent[i][1],t,h] - 2*(data.r[i]*P[i,t,h] + data.x[i]*Q[i,t,h]) + (data.r[i]^2 + data.x[i]^2)*l[i,t,h]
            
            [i in data.L, t in 1:data.T, h in 1:data.H],
            l[i,t,h]*v[data.parent[i][1],t,h] >= P[i,t,h]^2 + Q[i,t,h]^2
            
            [i in data.L, t in 1:data.T, h in 1:data.H],
            P[i,t,h]^2 + Q[i,t,h]^2 <= data.A[i]^2
            
            [i in data.L, t in 1:data.T, h in 1:data.H],
            (P[i,t,h]-data.r[i]*l[i,t,h])^2 + (Q[i,t,h]-data.x[i]*l[i,t,h])^2 <= data.A[i]^2   
        end
    )

    @expression(model, transportation_cost, sum(data.Ct[m,t]*z[m,t] for m in 1:data.M, t in 2:data.T))
    @expression(model, power_gen_cost, sum(data.alpha[t]*data.Cp[i,t,h]*pg[i,t,h] for h in 1:data.H, t in 1:data.T, i in data.B))
    @expression(model, VoLL, sum(data.alpha[t]*data.C_VoLL[i]*(1-delta[i,t,h])*pd[t][i,h] for h in 1:data.H, t in 1:data.T, i in data.B))
    @expression(model, genset_fuel_cost, sum(data.alpha[t]*data.Ca[t]*p_a[m,i,t,h] for h in 1:data.H, m in 1:data.M, t in 1:data.T, i in data.B))

    @objective(model, Min, transportation_cost + power_gen_cost + VoLL + genset_fuel_cost)

    optimize!(model)

    #-----------------------------CONVERTING y TO SPACE OF w---------------------------
    status = termination_status(model)
    solveTime = solve_time(model)

    objective_val = []
    for soln_index in 1:result_count(model)
        push!(objective_val, objective_value(model; result=soln_index) )
    end

    node_specific_w_value = []
    num_columns = result_count(model)
    for soln_index in 1:num_columns
        push!(node_specific_w_value, Dict{Int64, Matrix{Float64}}()) # node->column
        for t in 1:data.T
            for node in nodes_in_stage(tree, t+1)
                node_specific_w_value[soln_index][node] = zeros(data.M, length(data.B))
                for m in 1:data.M
                    for i in data.B
                        node_specific_w_value[soln_index][node][m,i] = round(value(y[m,i,t]; result=soln_index))
                    end
                end
            end
        end
    end
    
    
    return ExpectedValueSolution(status, num_columns, node_specific_w_value, objective_val, solveTime)

end


# Assuming transition probabilities of children nodes of a node `n` are all same
function expected_power_demand(data::ModelDataNodalStochastic; scenario_tree)
    
    approximate_pd = Matrix{Float64}[]
    approximate_qd = Matrix{Float64}[]
    num_of_buses = length(data.B)
    for stage in 2:data.T+1
        pd = zeros(num_of_buses, data.H)
        qd = zeros(num_of_buses, data.H)
        number_of_nodes = 0
        for node in nodes_in_stage(scenario_tree, stage)
            number_of_nodes += 1
            pd .+= data.pd[node]
            qd .+= data.qd[node] 
        end
        pd ./= number_of_nodes
        qd ./= number_of_nodes
        push!(approximate_pd, round.(pd; digits=5))
        push!(approximate_qd, round.(qd; digits=5))
    end

    return approximate_pd, approximate_qd # correspond to time periods and not stages
end

function index_of_buses_with_max_load(data::ModelDataNodalStochastic; tree::Tree)
    apd = expected_power_demand(data; scenario_tree=tree) 
    overall_load = sum(apd[1][1], dims=2)
    for t in 2:data.T
        overall_load .+= sum(apd[1][t], dims=2)
    end
    
    return sortperm(vec(overall_load), rev=true)
end