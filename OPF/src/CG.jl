function get_cost_for_fixed_w(data::ModelDataNodalStochastic;
    tree::Tree, 
    nfix::Int64,
    w_values::Matrix{Float64}
    )

    @assert nfix > 1

    for m in 1:data.M
        located_at_buses = 0
        for i in data.B
            located_at_buses += w_values[m,i]
        end
        @assert located_at_buses == 1 "Invalid column!"
    end

    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "TimeLimit", 7200.0)

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
        p_a[m,i,h] <= data.pa[m]*w_values[m,i]
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

    @objective(model, Min, power_gen_cost + VoLL + genset_fuel_cost)

    optimize!(model)

    #-----------------------------RETRIEVING SOLUTION---------------------------
    solution = SolutionInfo() 
    solution.status = termination_status(model)
    solution.solvetime = solve_time(model)
    if solution.status == MOI.INFEASIBLE || solution.status == MOI.DUAL_INFEASIBLE || solution.status == MOI.INFEASIBLE_OR_UNBOUNDED
        solution.objective = "NA"
        solution.bound = "NA"
        solution.gap = "NA"
        solution.sol = "NA"
    else
        solution.objective = objective_value(model)
        solution.decisions = (value.(pg), value.(qg), value.(delta), value.(v), value.(l), value.(p_a), value.(q_a), value.(P), value.(Q))
    end
    #---------------------------------------------------------------------------

    solution.num_variables = num_variables(model)
    solution.num_constraints = num_constraints(model; count_variable_in_set_constraints = true)
    solution.num_structural_constraints = num_constraints(model; count_variable_in_set_constraints = false)

    return solution
end

function solve_CG(data::ModelDataNodalStochastic; 
    tree::Tree,
    iterlimit::Int64 = 1000,
    timelimit::Float64 = 7200.0,
    relative_gap::Float64 = 1e-6,
    sp_gap = 1e-4,
    do_column_sharing::Bool = false)

    CG_solution = CGSolutionInfo()
    CG_solution.SP_solvetime_per_iteration = Float64[]
    CG_solution.CS_solvetime_per_iteration = Float64[]
    CG_solution.RMP_solvetime_per_iteration = Float64[]
    
    N = maximum(nodes(tree))
    total_solution_time = 0.0
    total_solution_time_perfect_parallelization = 0.0

    if do_column_sharing 
        sibling_nodes_pairs = pairs_of_sibling_nodes(tree)
    end

    total_buses = length(data.B)

    #------------------------Initial columns and corresponding cost---------------------------
    # Generating initial columns
    evp_solution = deterministic_model(data; tree = tree, timelimit = 300, gap = 1e-4, desired_num_sols=50)
    total_solution_time += evp_solution.solvetime
    total_solution_time_perfect_parallelization += evp_solution.solvetime

    # NC contains number of columns for each subproblem
    NC = ones(Int, N)
    # Dictionary that will contain columns priced out in each subproblem.
    columns = Dict()

    # no columns for the root node, not a subproblem
    NC[1] = 0
    columns[1] = []
    
    idx_bus = index_of_buses_with_max_load(data; tree=tree)
    for n in 2:N columns[n] = [] end
    for soln_index in 0:evp_solution.num_of_columns
        for n in 2:N
            if soln_index == 0 
                push!(columns[n], zeros(data.M, total_buses)) 
                for m in 1:data.M
                    columns[n][1][m,idx_bus[m]] = 1
                end
            else
                if evp_solution.column[soln_index][n] in columns[n]
                    continue        # Don't add column if it's already added
                else
                    NC[n] += 1
                    push!(columns[n], evp_solution.column[soln_index][n])
                end
            end
        end
    end
    
    column_cost = Dict()
    column_cost[1] = []

    pair = [] # to store the pair of (n, column) for parallelization
    index_tracker = 0
    for n in 2:N
        for col_index in 1:NC[n]
            index_tracker += 1
            push!(pair, (n, col_index))
        end
    end

    ICG_total_elapsed_time = @elapsed begin
        z = pmap(pair->get_cost_for_fixed_w(data; tree= tree, nfix = pair[1], w_values = columns[pair[1]][pair[2]]), pair; retry_delays = zeros(3))
    end

    time_holder = Float64[]
    for n in 2:N column_cost[n] = [] end
    for idx in 1:index_tracker
        push!(column_cost[pair[idx][1]], z[idx].objective)
        push!(time_holder, z[idx].solvetime)
    end
    total_solution_time += ICG_total_elapsed_time
    total_solution_time_perfect_parallelization += maximum(time_holder)
    CG_solution.initial_col_gen_time = total_solution_time
    #----------------------------------------------------------------------------

    #------------------initializing relative gap tracker------------------
    relGapTracker = []
    #---------------------------------------------------------------------------

    #----------------------------INITIALIZING BOUNDS----------------------------
    LB = -1e+10
    UB = +1e+10
    #---------------------------------------------------------------------------

    iter = 0
    cs_iter = 0 # number of iterations where CS is applied
    rmp_solution = RMPSolutionInfo()
    while true
        
        iter += 1

        #-----------------------SOLVING RELAXED RMP-----------------------------
        rmp_solution = solve_RMP(data, column_cost, columns, NC; tree = tree, relaxed = true)
        total_solution_time += rmp_solution.solvetime
        total_solution_time_perfect_parallelization += rmp_solution.solvetime
        push!(CG_solution.RMP_solvetime_per_iteration, rmp_solution.solvetime)

        @assert (rmp_solution.objective <= (UB + 1e-2) ) "UB goes up from $(UB) to $(rmp_solution.objective) !"

        UB = rmp_solution.objective
        #-----------------------------------------------------------------------

        #-----------------------------SOLVING SP--------------------------------
        sp_solution = SPSolutionInfo[]
        local_LB = 0.0
        time_holder = Float64[]

        SP_total_elapsed_time = @elapsed begin 
            z = pmap(n->solve_SP(data, n, rmp_solution; tree = tree, gap_tol = sp_gap), 2:N; retry_delays = zeros(3))
        end
        for n in 1:N
            if n == 1
                push!(sp_solution, SPSolutionInfo())
            else
                push!(sp_solution, z[n-1])
                push!(time_holder, sp_solution[n].solvetime)
                local_LB += sp_solution[n].bound
            end
        end
        total_solution_time += SP_total_elapsed_time
        total_solution_time_perfect_parallelization += maximum(time_holder)
        push!(CG_solution.SP_solvetime_per_iteration, SP_total_elapsed_time)

        #-------Updating LB---------------#
        LB = max(LB, (UB + local_LB ) ) 
        #---------------------------------#

        #--------Calculating relative gap------#
        relGap = abs(UB-LB)/((1e-10)+abs(UB))
        push!(relGapTracker, relGap)
        #--------------------------------------#
        
        #---------------------CHECKING TERMINATION CRITERIA------------------------------
        if relGap <= relative_gap || iter == iterlimit || total_solution_time >= timelimit
            break 
        end
        #-----------------------------------------------------------------------

        #----------------Update columns to be added, corresponding cost and NC---------------#
        for n in 2:N
            for col in sp_solution[n].columns_info
                if col.add_column
                    NC[n] += 1
                    push!(columns[n], col.column)
                    push!(column_cost[n], col.cost)
                end
            end
        end
        #------------------------------------------------------------------------------------#
        
        #-----------------Column sharing (CS) procedure--------------------#
        if do_column_sharing
            sibling_pairs_with_col_info = pairs_of_sibling_nodes(tree, sibling_nodes_pairs, sp_solution, columns)
            CS_total_time_elapsed = @elapsed begin
                CS_sp_solution = pmap(field->column_sharing(data,field[1], field[2], field[3];
                    tree = tree, sp_core_solution = sp_solution),
                    sibling_pairs_with_col_info; retry_delays = zeros(3))
            end

            # Updating columns from CS procedure
            time_holder = Float64[]
            for (index, pair) in enumerate(sibling_pairs_with_col_info)
                target_node = pair[2]
                if CS_sp_solution[index].columns_info.add_column
                    NC[target_node] += 1
                    push!(columns[target_node], CS_sp_solution[index].columns_info.column)
                    push!(column_cost[target_node], CS_sp_solution[index].columns_info.cost)
                end
                push!(time_holder, CS_sp_solution[index].solvetime)
            end

            total_solution_time += CS_total_time_elapsed
            total_solution_time_perfect_parallelization += maximum(time_holder; init=0)
            push!(CG_solution.CS_solvetime_per_iteration, CS_total_time_elapsed)

            cs_iter += 1
        end
        #-------------------------------------------------------------#
    
    
    end
    
    # Check if the CG solution is integral
    is_y_integral = all(is_solution_binary.(rmp_solution.y_val))
    is_z_integral = all(is_solution_binary.(rmp_solution.z_val))
    is_rho_integral = all(is_solution_binary.(rmp_solution.rho_val))



    # if not integral solution, then solve MIP in the end
    if !is_y_integral || !is_rho_integral || !is_z_integral
        CG_solution.solution_integral = false
        CG_solution.MIP_RMP_solution = solve_RMP(data, column_cost, columns, NC; tree = tree, relaxed = false)
        total_solution_time += CG_solution.MIP_RMP_solution.solvetime
        total_solution_time_perfect_parallelization += CG_solution.MIP_RMP_solution.solvetime
    end

    CG_solution.UB = UB
    CG_solution.LB = LB
    CG_solution.relaxed_RMP_solution = rmp_solution
    CG_solution.solvetime = total_solution_time
    CG_solution.solvetime_perfect_parallelization = total_solution_time_perfect_parallelization
    CG_solution.iterations = iter
    CG_solution.relative_gap_each_iter = relGapTracker
    CG_solution.final_relative_gap = relGapTracker[end]
    CG_solution.column_generated = columns
    CG_solution.columns_cost = column_cost
    CG_solution.total_columns_generated = sum(NC)
    CG_solution.CS_iterations = cs_iter

    return CG_solution
end


