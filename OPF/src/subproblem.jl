function solve_SP(data::ModelDataNodalStochastic,
    nfix, 
    RMP_solution;
    tree,
    gap_tol=1e-4)

    @assert nfix > 1

    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "TimeLimit", 7200.0)
    set_optimizer_attribute(model, "MIPGap", gap_tol)
    set_optimizer_attribute(model, "PoolSolutions", 5)
    set_optimizer_attribute(model, "PoolSearchMode", 2)

    #-----------------------------CONSTRUCT MODEL-------------------------------
    @variables(
        model,
        begin
            w[m in 1:data.M, i in data.B], Bin

            0 <= pg[i in data.B, h in 1:data.H] <= data.pg_ub[i]
            0 <= qg[i in data.B, h in 1:data.H] <= data.qg_ub[i]
            0 <= delta[i in data.B, h in 1:data.H] <= 1
            data.v_lb[i] <= v[i in data.B, h in 1:data.H] <= data.v_ub[i]
            l[i in data.L, h in 1:data.H] >= 0
            p_a[m in 1:data.M, i in data.B, h in 1:data.H] >= 0
            
            q_a[m in 1:data.M, i in data.B, h in 1:data.H]
            P[i in data.L, h in 1:data.H]
            Q[i in data.L, h in 1:data.H]
        end
    )

    @constraint(
        model,
        [m in 1:data.M], sum(w[m, i] for i in data.B) == 1
    )

    # Power generation constraint
    @constraint(
        model,
        [m in 1:data.M, i in data.B, h in 1:data.H],
        p_a[m,i,h] <= data.pa[m]*w[m,i]
    )
    @constraints(model, begin
        [m in 1:data.M, i in data.B, h in 1:data.H], q_a[m,i,h] <= p_a[m,i,h] 
        [m in 1:data.M, i in data.B, h in 1:data.H], q_a[m,i,h] >= -p_a[m,i,h]
    end)

    # Power flow equations
    @constraints(
        model,
        begin
            [i in data.B, h in 1:data.H],
                pg[i,h] - data.pd[nfix][i,h]*delta[i,h] + sum(p_a[m,i,h] for m in 1:data.M) ==
                sum(P[j,h] for j in data.children[i]) - sum((P[i,h]-data.r[i]*l[i,h]) for _ in data.parent[i]) + data.g[i]*v[i,h] 
            
            [i in data.B, h in 1:data.H],
            qg[i,h] - data.qd[nfix][i,h]*delta[i,h] + sum(q_a[m,i,h] for m in 1:data.M) ==
            sum(Q[j,h] for j in data.children[i]) - sum((Q[i,h]-data.x[i]*l[i,h]) for _ in data.parent[i]) + data.b[i]*v[i,h] 
        end
    )

    @constraints(
        model, 
        begin
            [i in data.L, h in 1:data.H],
            v[i,h] == v[data.parent[i][1],h] - 2*(data.r[i]*P[i,h] + data.x[i]*Q[i,h]) + (data.r[i]^2 + data.x[i]^2)*l[i,h]
            
            [i in data.L, h in 1:data.H],
            l[i,h]*v[data.parent[i][1],h] >= P[i,h]^2 + Q[i,h]^2
            
            [i in data.L, h in 1:data.H],
            P[i,h]^2 + Q[i,h]^2 <= data.A[i]^2
            
            [i in data.L, h in 1:data.H],
            (P[i,h]-data.r[i]*l[i,h])^2 + (Q[i,h]-data.x[i]*l[i,h])^2 <= data.A[i]^2   
        end
    )

    #define objective terms
    @expression(model, power_gen_cost, data.prob[nfix]*sum(data.alpha[stage(tree,nfix)[1]]*data.Cp[i,stage(tree,nfix)[1],h]*pg[i,h] for h in 1:data.H, i in data.B))
    @expression(model, VoLL, data.prob[nfix]*sum(data.alpha[stage(tree,nfix)[1]]*data.C_VoLL[i]*(1-delta[i,h])*data.pd[nfix][i,h] for h in 1:data.H, i in data.B))
    @expression(model, genset_fuel_cost, data.prob[nfix]*sum(data.alpha[stage(tree,nfix)[1]]*data.Ca[stage(tree,nfix)[1]]*p_a[m,i,h] for h in 1:data.H, m in 1:data.M, i in data.B))
    
    if nfix in leaves(tree)[1]
        @expression(model, dual_terms, sum(RMP_solution.gamma1[m,i,nfix]*w[m,i] for m in 1:data.M, i in data.B) + RMP_solution.mu[nfix])
    else
        @expression(model, dual_terms, sum(RMP_solution.gamma1[m,i,nfix]*w[m,i] for m in 1:data.M, i in data.B) -
                                       sum(RMP_solution.gamma2[m,i,nfix]*w[m,i] for m in 1:data.M, i in data.B) +
                                       sum(RMP_solution.gamma3[m,i,nfix]*w[m,i] for m in 1:data.M, i in data.B) +
                                       RMP_solution.mu[nfix])
    end

    @objective(model, Min, power_gen_cost + VoLL + genset_fuel_cost - dual_terms)
    #---------------------------------------------------------------------------
    
    optimize!(model)

    SP_solution = SPSolutionInfo()
    #-----------------------------RETRIEVING STATS---------------------------
    SP_solution.status = termination_status(model)
    SP_solution.solvetime = solve_time(model)
    SP_solution.bound = objective_bound(model)
    #---------------------------------------------------------------------------

    #---------------------RETURNING COLUMNS TO BE ADDED (IF APPLICABLE)---------
    SP_solution.columns_info = []
    for soln_index in 1:result_count(model)
        push!(SP_solution.columns_info, SPColumnInfo())
        if objective_value(model; result = soln_index) < -1e-5
            SP_solution.columns_info[soln_index].add_column = true
            SP_solution.columns_info[soln_index].column = round.(value.(w; result = soln_index))
            SP_solution.columns_info[soln_index].cost = objective_value(model; result = soln_index) + value(dual_terms; result = soln_index)
        end
    end
    #---------------------------------------------------------------------------

    return SP_solution
end

function column_sharing(data::ModelDataNodalStochastic, 
    source_node::Int64, 
    target_node::Int64,
    column_number::Int64; 
    tree::Tree,
    sp_core_solution::Vector{SPSolutionInfo},
    timelimit::Float64 = 7200.0)

    @assert (source_node > 1 && target_node > 1) "Invalid column sharing!"

    new_sp_solution = SPSolutionInfo()
    new_sp_solution.solvetime = 0.0
    new_sp_solution.columns_info = SPColumnInfo()

    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "TimeLimit", timelimit)

    @variables(
        model,
        begin
            0 <= pg[i in data.B, h in 1:data.H] <= data.pg_ub[i]
            0 <= qg[i in data.B, h in 1:data.H] <= data.qg_ub[i]
            0 <= delta[i in data.B, h in 1:data.H] <= 1
            data.v_lb[i] <= v[i in data.B, h in 1:data.H] <= data.v_ub[i]
            l[i in data.L, h in 1:data.H] >= 0
            p_a[m in 1:data.M, i in data.B, h in 1:data.H] >= 0
            
            q_a[m in 1:data.M, i in data.B, h in 1:data.H]
            P[i in data.L, h in 1:data.H]
            Q[i in data.L, h in 1:data.H]
        end
    )

    # Power generation constraint
    @constraint(
        model,
        [m in 1:data.M, i in data.B, h in 1:data.H],
        p_a[m,i,h] <= data.pa[m]*(sp_core_solution[source_node].columns_info[column_number].column[m,i])
    )
    @constraints(model, begin
        [m in 1:data.M, i in data.B, h in 1:data.H], q_a[m,i,h] <= p_a[m,i,h] 
        [m in 1:data.M, i in data.B, h in 1:data.H], q_a[m,i,h] >= -p_a[m,i,h]
    end)

    # Power flow equations
    @constraints(
        model,
        begin
            [i in data.B, h in 1:data.H],
                pg[i,h] - data.pd[target_node][i,h]*delta[i,h] + sum(p_a[m,i,h] for m in 1:data.M) ==
                sum(P[j,h] for j in data.children[i]) - sum((P[i,h]-data.r[i]*l[i,h]) for _ in data.parent[i]) + data.g[i]*v[i,h] 
            
            [i in data.B, h in 1:data.H],
            qg[i,h] - data.qd[target_node][i,h]*delta[i,h] + sum(q_a[m,i,h] for m in 1:data.M) ==
            sum(Q[j,h] for j in data.children[i]) - sum((Q[i,h]-data.x[i]*l[i,h]) for _ in data.parent[i]) + data.b[i]*v[i,h] 
        end
    )

    @constraints(
        model, 
        begin
            [i in data.L, h in 1:data.H],
            v[i,h] == v[data.parent[i][1],h] - 2*(data.r[i]*P[i,h] + data.x[i]*Q[i,h]) + (data.r[i]^2 + data.x[i]^2)*l[i,h]
            
            [i in data.L, h in 1:data.H],
            l[i,h]*v[data.parent[i][1],h] >= P[i,h]^2 + Q[i,h]^2
            
            [i in data.L, h in 1:data.H],
            P[i,h]^2 + Q[i,h]^2 <= data.A[i]^2
            
            [i in data.L, h in 1:data.H],
            (P[i,h]-data.r[i]*l[i,h])^2 + (Q[i,h]-data.x[i]*l[i,h])^2 <= data.A[i]^2   
        end
    )

    #define objective terms
    @expression(model, power_gen_cost, data.prob[target_node]*sum(data.alpha[stage(tree,target_node)[1]]*data.Cp[i,stage(tree,target_node)[1],h]*pg[i,h] for h in 1:data.H, i in data.B))
    @expression(model, VoLL, data.prob[target_node]*sum(data.alpha[stage(tree,target_node)[1]]*data.C_VoLL[i]*(1-delta[i,h])*data.pd[target_node][i,h] for h in 1:data.H, i in data.B))
    @expression(model, genset_fuel_cost, data.prob[target_node]*sum(data.alpha[stage(tree,target_node)[1]]*data.Ca[stage(tree,target_node)[1]]*p_a[m,i,h] for h in 1:data.H, m in 1:data.M, i in data.B))

    @objective(model, Min, power_gen_cost + VoLL + genset_fuel_cost) 
    
    optimize!(model)

    #-----------------------------RETRIEVING SOLUTION---------------------------
    new_sp_solution.status = termination_status(model)
    new_sp_solution.solvetime = solve_time(model)

    if result_count(model) >= 1
        new_sp_solution.columns_info.add_column = true
        new_sp_solution.columns_info.column = sp_core_solution[source_node].columns_info[column_number].column
        new_sp_solution.columns_info.cost = objective_value(model) 
    end
    #---------------------------------------------------------------------------


    return new_sp_solution
end