function get_cost_for_fixed_w(data::ModelDataNodalStochastic;
    tree::Tree, 
    nfix::Int64,
    w_values::Vector{Float64})
    @assert nfix > 1

    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "NonConvex", 2)
    set_optimizer_attribute(model, "TimeLimit", 7200.0)

    # define variables
    @variables(
        model,
        begin
            data.d_min[nfix][j] <= d[j in 1:data.J] <= data.d_max[nfix][j]
            0 <= c[j in 1:data.J] <= 1
            F[i in 1:data.I, j in 1:data.J] >= 0
        end
    )

    # Feasible production
    @constraint(
        model,
        [i in 1:data.I],
        sum(F[i, j] for j in 1:data.J) <= data.C[i]*round(w_values[i])
    )

    # demand satisfcation
    @constraint(
        model,
        [j in 1:data.J],
        sum(F[i, j] for i in 1:data.I) == d[j]
    )

    # blending
    @constraint(
        model,
        [j in 1:data.J],
        sum(data.λ[i] * F[i, j] for i in 1:data.I) == c[j] * d[j]
    )

    @expression(model, expr_b, data.prob[nfix] * sum(data.b[i,stage(tree, nfix)[1]] * F[i,j] for j in 1:data.J, i in 1:data.I) )
    @expression(model, expr_r, data.prob[nfix] * sum(data.r[i,j]*F[i,j] for i in 1:data.I, j in 1:data.J) )
    @expression(model, expr_F, data.prob[nfix] * sum((data.price_slope[j, stage(tree, nfix)[1]] * c[j] + data.price_intercept[j, stage(tree, nfix)[1]]) * F[i,j] for i in 1:data.I, j in 1:data.J) )
    
    # define objective
    @objective(model, Min, expr_r + expr_b - expr_F) 
    
    optimize!(model)

    #-----------------------------RETRIEVING SOLUTION---------------------------
    status = termination_status(model)
    solveTime = solve_time(model)

    if status == MOI.INFEASIBLE || status == MOI.DUAL_INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
        objective = "NA"
        bound = "NA"
        gap = "NA"
        sol = "NA"
    else
        objective = objective_value(model)
        bound = objective_bound(model)
        gap = abs(objective - bound) / ((1e-10) + abs(objective))
        sol = (value.(d), value.(c), value.(F))
    end
    #---------------------------------------------------------------------------

    num_of_vars = num_variables(model)
    num_of_constraints = num_constraints(model; count_variable_in_set_constraints = true)
    num_of_structural_constraints = num_constraints(model; count_variable_in_set_constraints = false)

    return SolutionInfo(status, objective, bound, solveTime, gap, sol, num_of_vars, num_of_constraints, num_of_structural_constraints)
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

    #------------------------Initial columns and corresponding cost---------------------------
    # NC contains number of columns for each subproblem (including root node; although, it is not a subproblem)
    # We start with one column per subproblem 
    NC = ones(Int, N)

    # Dictionary that will contain columns priced out in each subproblem. Again, for completeness I have a key corrsponding 
    # to root node; although, it is not a subproblem here.
    columns = Dict()
    for n in 1:N 
        if n == 1 
            NC[n] = 0 # Set number of columns for root node to zero
            columns[n] = Float64[]
        else 
            columns[n] = reshape(ones(data.I), 1, data.I)
        end
    end

    
    column_cost = Dict()
    column_cost[1] = Float64[] # Set column cost to be an empty vector for root node
    
    ICG_total_elapsed_time = @elapsed begin
        z = pmap(n->get_cost_for_fixed_w(data; tree= tree, nfix = n, w_values = columns[n][1,:]), 2:N; retry_delays = zeros(3))
    end
    time_holder = Float64[]
    for n in 2:N
        column_cost[n] = [z[n-1].objective]
        push!(time_holder, z[n-1].solvetime)
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
        for n in 1:N
            if sp_solution[n].add_column
                NC[n] += 1
                columns[n] = vcat(columns[n], reshape(sp_solution[n].column, 1, data.I))
                push!(column_cost[n], sp_solution[n].column_cost)
            end
        end
        #------------------------------------------------------------------------------------#
        
        #-----------------Column sharing (CS) procedure--------------------#
        if do_column_sharing
            sibling_nodes_pairs = pairs_of_sibling_nodes(tree)
            CS_total_time_elapsed = @elapsed begin
                CS_sp_solution = pmap(nodes->column_sharing(data,nodes[1], nodes[2];
                    tree = tree, sp_core_solution = sp_solution),
                    sibling_nodes_pairs; retry_delays = zeros(3))
            end

            # Updating columns from CS procedure
            time_holder = Float64[]
            for (index, pair) in enumerate(sibling_nodes_pairs)
                target_node = pair[2]
                if CS_sp_solution[index].add_column
                    NC[target_node] += 1
                    columns[target_node] = vcat(columns[target_node], reshape(CS_sp_solution[index].column, 1, data.I))
                    push!(column_cost[target_node], CS_sp_solution[index].column_cost)
                end
                push!(time_holder, CS_sp_solution[index].solvetime)
            end

            total_solution_time += CS_total_time_elapsed
            total_solution_time_perfect_parallelization += maximum(time_holder)
            push!(CG_solution.CS_solvetime_per_iteration, CS_total_time_elapsed)
            
            cs_iter += 1
        end
        #-------------------------------------------------------------#
    end
    
    # Check if the CG solution is integral, i.e, rho and x 0/1 or not?
    is_x_integral = all(is_solution_binary.(rmp_solution.x_val))
    is_rho_integral = all(is_solution_binary.(rmp_solution.rho_val))

    # if not integral solution, then solve MIP in the end
    if !is_x_integral || !is_rho_integral
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


