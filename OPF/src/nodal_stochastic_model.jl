function solve_stoch_model_FS(data::ModelDataNodalStochastic;
    tree::Tree,
    timelimit,
    gap
)

    N = maximum(nodes(tree))

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", timelimit)
    set_optimizer_attribute(model, "MIPGap", gap)

    @variables(
        model,
        begin
            # binary variables
            y[m in 1:data.M, i in data.B, n in 1:(minimum(leaves(tree)[1])-1)], Bin
            z[m in 1:data.M, n in 2:(minimum(leaves(tree)[1])-1)], Bin

            # (semi-)bounded continuous variables
            0 <= pg[i in data.B, n in 2:N, h in 1:data.H] <= data.pg_ub[i]
            0 <= qg[i in data.B, n in 2:N, h in 1:data.H] <= data.qg_ub[i]
            0 <= delta[i in data.B, n in 2:N, h in 1:data.H] <= 1
            data.v_lb[i] <= v[i in data.B, n in 2:N, h in 1:data.H] <= data.v_ub[i]
            l[i in data.L, n in 2:N, h in 1:data.H] >= 0
            p_a[m in 1:data.M, i in data.B, n in 2:N, h in 1:data.H] >= 0

            # unrestricted
            q_a[m in 1:data.M, i in data.B, n in 2:N, h in 1:data.H]
            P[i in data.L, n in 2:N, h in 1:data.H]
            Q[i in data.L, n in 2:N, h in 1:data.H]
        end
    )

    @constraint(
        model,
        [m in 1:data.M, n in 1:(minimum(leaves(tree)[1])-1)], sum(y[m, i, n] for i in data.B) == 1
    )

    @constraints(
        model,
        begin
            [m in 1:data.M, i in data.B, n in 2:(minimum(leaves(tree)[1])-1)], y[m,i,n]-y[m,i,tree.parent[n]] <= z[m,n]
            [m in 1:data.M, i in data.B, n in 2:(minimum(leaves(tree)[1])-1)], z[m,n] <= 2 - y[m,i,n] - y[m,i,tree.parent[n]] 
        end
    )

    # Power generation constraint
    @constraint(
        model,
        [m in 1:data.M, i in data.B, n in 2:N, h in 1:data.H],
        p_a[m,i,n,h] <= data.pa[m]*y[m,i,tree.parent[n]]
    )
    @constraints(model, begin
        [m in 1:data.M, i in data.B, n in 2:N, h in 1:data.H], q_a[m,i,n,h] <= p_a[m,i,n,h] 
        [m in 1:data.M, i in data.B, n in 2:N, h in 1:data.H], q_a[m,i,n,h] >= -p_a[m,i,n,h]
    end)

    # Power flow equations
    @constraints(
        model,
        begin
            [i in data.B, n in 2:N, h in 1:data.H],
                pg[i,n,h] - data.pd[n][i,h]*delta[i,n,h] + sum(p_a[m,i,n,h] for m in 1:data.M) ==
                sum(P[j,n,h] for j in data.children[i]) - sum((P[i,n,h]-data.r[i]*l[i,n,h]) for _ in data.parent[i]) + data.g[i]*v[i,n,h] 
            
            [i in data.B, n in 2:N, h in 1:data.H],
            qg[i,n,h] - data.qd[n][i,h]*delta[i,n,h] + sum(q_a[m,i,n,h] for m in 1:data.M) ==
            sum(Q[j,n,h] for j in data.children[i]) - sum((Q[i,n,h]-data.x[i]*l[i,n,h]) for _ in data.parent[i]) + data.b[i]*v[i,n,h] 
        end
    )

    @constraints(
        model, 
        begin
            [i in data.L, n in 2:N, h in 1:data.H],
            v[i,n,h] == v[data.parent[i][1],n,h] - 2*(data.r[i]*P[i,n,h] + data.x[i]*Q[i,n,h]) + (data.r[i]^2 + data.x[i]^2)*l[i,n,h]
            
            [i in data.L, n in 2:N, h in 1:data.H],
            l[i,n,h]*v[data.parent[i][1],n,h] >= P[i,n,h]^2 + Q[i,n,h]^2
            
            [i in data.L, n in 2:N, h in 1:data.H],
            P[i,n,h]^2 + Q[i,n,h]^2 <= data.A[i]^2
            
            [i in data.L, n in 2:N, h in 1:data.H],
            (P[i,n,h]-data.r[i]*l[i,n,h])^2 + (Q[i,n,h]-data.x[i]*l[i,n,h])^2 <= data.A[i]^2
        end
    )

    #define objective terms
    @expression(model, transportation_cost, sum(data.prob[n]*data.Ct[m,stage(tree,n)[1]+1]*z[m,n] for m in 1:data.M, n in 2:(minimum(leaves(tree)[1])-1)))
    @expression(model, power_gen_cost, sum(data.prob[n]*data.alpha[stage(tree,n)[1]]*data.Cp[i,stage(tree,n)[1],h]*pg[i,n,h] for h in 1:data.H, i in data.B, n in 2:N))
    @expression(model, VoLL, sum(data.prob[n]*data.alpha[stage(tree,n)[1]]*data.C_VoLL[i]*(1-delta[i,n,h])*data.pd[n][i,h] for h in 1:data.H, i in data.B, n in 2:N))
    @expression(model, genset_fuel_cost, sum(data.prob[n]*data.alpha[stage(tree,n)[1]]*data.Ca[stage(tree,n)[1]]*p_a[m,i,n,h] for h in 1:data.H, m in 1:data.M, i in data.B, n in 2:N))

    @objective(model, Min, transportation_cost + power_gen_cost + VoLL + genset_fuel_cost)

    optimize!(model)

    #-----------------------------RETRIEVING SOLUTION---------------------------
    status = termination_status(model)
    solveTime = solve_time(model)

    if result_count(model) == 0
        objective = "NA"
        bound = "NA"
        gap = "NA"
        sol = "NA"
    else
        objective = objective_value(model)
        bound = objective_bound(model)
        gap = abs(objective - bound) / ((1e-10) + abs(objective))
        sol = (value.(y), value.(z), value.(pg), value.(qg), value.(delta), value.(v), value.(l), value.(p_a), value.(q_a), value.(P), value.(Q))
    end

    num_of_vars = num_variables(model)
    num_of_constraints = num_constraints(model; count_variable_in_set_constraints = true)
    num_of_structural_constraints = num_constraints(model; count_variable_in_set_constraints = false)
    #---------------------------------------------------------------------------

    return SolutionInfo(status, objective, bound, solveTime, gap, sol, num_of_vars, num_of_constraints, num_of_structural_constraints)
end
