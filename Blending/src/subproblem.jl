# Args: data - Model data
#       nfix - Node index (corresponds to the the subproblem for this node)
#       RMP_solution - Solution of the RMP
#       tree - Scenario tree
#       gap_tol - Gap tolerance
# Returns: SP_solution - Solution of the subproblem
function solve_SP(data::ModelDataNodalStochastic,
    nfix,
    RMP_solution;
    tree,
    gap_tol=1e-4)

    @assert nfix > 1 # There is no subproblem corresponding to the root node

    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "NonConvex", 2)
    set_optimizer_attribute(model, "TimeLimit", 7200.0)
    set_optimizer_attribute(model, "MIPGap", gap_tol)

    #-----------------------------CONSTRUCT MODEL-------------------------------
    @variables(
        model,
        begin
            w[i in 1:data.I], Bin
            0 <= c[j in 1:data.J] <= 1
            data.d_min[nfix][j] <= d[j in 1:data.J] <= data.d_max[nfix][j]
            F[i in 1:data.I, j in 1:data.J] >= 0
        end
    )

    # Feasible production
    @constraint(
        model,
        [i in 1:data.I],
        sum(F[i, j] for j in 1:data.J) <= data.C[i]*w[i]
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
    
    if nfix in leaves(tree)[1]
        @expression(model, dual_terms, sum(RMP_solution.gamma2[i,nfix]*w[i] for i in 1:data.I) + RMP_solution.mu[nfix])
    else
        @expression(model, dual_terms, sum((RMP_solution.gamma1[i,nfix]+RMP_solution.gamma2[i,nfix])*w[i] for i in 1:data.I) + RMP_solution.mu[nfix])
    end

    @objective(model, Min, expr_r + expr_b - expr_F - dual_terms)
    #---------------------------------------------------------------------------
    
    optimize!(model)

    SP_solution = SPSolutionInfo()
    #-----------------------------RETRIEVING STATUS---------------------------
    SP_solution.status = termination_status(model)
    SP_solution.solvetime = solve_time(model)
    #---------------------------------------------------------------------------

    #-----------------------------RETRIEVING SOLUTION---------------------------
    SP_solution.objective = objective_value(model)
    SP_solution.bound = objective_bound(model)
    SP_solution.gap = abs(SP_solution.objective - SP_solution.bound) / ((1e-10) + abs(SP_solution.objective))
    #---------------------------------------------------------------------------

    #---------------------RETURNING COLUMNS TO BE ADDED (IF APPLICABLE)---------
    if SP_solution.objective < -1e-5
        SP_solution.add_column = true
        SP_solution.column = round.(value.(w))
        SP_solution.column_cost = SP_solution.objective + value(dual_terms)
    end
    #---------------------------------------------------------------------------

    return SP_solution
end

# Args: data - Model data
#       source_node - Source node index
#       target_node - Target node index
#       tree - Scenario tree
#       sp_core_solution - vector of original subproblem solutions
#       timelimit - Time limit for solving the subproblem
function column_sharing(data::ModelDataNodalStochastic, 
    source_node::Int64, 
    target_node::Int64; 
    tree::Tree,
    sp_core_solution::Vector{SPSolutionInfo},
    timelimit::Float64 = 7200.0)

    @assert (source_node > 1 && target_node > 1) "Invalid column sharing!"

    new_sp_solution = SPSolutionInfo() # Populate this with the new column later on
    new_sp_solution.solvetime = 0.0

    # Proceed only if the target and source nodes have different decisions
    # Else exit, since sharing columns will result in the same decisions
    if !sp_core_solution[source_node].add_column
        return new_sp_solution
    elseif sp_core_solution[target_node].add_column
        if sp_core_solution[target_node].column == sp_core_solution[source_node].column
            return new_sp_solution
        end
    end


    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "NonConvex", 2)
    set_optimizer_attribute(model, "TimeLimit", timelimit)

    # define variables
    @variables(
        model,
        begin
            data.d_min[target_node][j] <= d[j in 1:data.J] <= data.d_max[target_node][j]
            0 <= c[j in 1:data.J] <= 1
            F[i in 1:data.I, j in 1:data.J] >= 0
        end
    )

    # Feasible production
    @constraint(
        model,
        [i in 1:data.I],
        sum(F[i, j] for j in 1:data.J) <= data.C[i]*(sp_core_solution[source_node].column[i])
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

    @expression(model, expr_b, data.prob[target_node] * sum(data.b[i,stage(tree, target_node)[1]] * F[i,j] for j in 1:data.J, i in 1:data.I) )
    @expression(model, expr_r, data.prob[target_node] * sum(data.r[i,j]*F[i,j] for i in 1:data.I, j in 1:data.J) )
    @expression(model, expr_F, data.prob[target_node] * sum( (data.price_slope[j, stage(tree, target_node)[1]] * c[j] + data.price_intercept[j, stage(tree, target_node)[1]]) * F[i,j] for i in 1:data.I, j in 1:data.J) )
    
    # define objective
    @objective(model, Min, expr_r + expr_b - expr_F) 
    
    optimize!(model)

    #-----------------------------RETRIEVING SOLUTION---------------------------
    new_sp_solution.status = termination_status(model)
    new_sp_solution.solvetime = solve_time(model)

    if result_count(model) >= 1
        new_sp_solution.objective = objective_value(model)
        new_sp_solution.bound = objective_bound(model)
        new_sp_solution.gap = abs(new_sp_solution.objective - new_sp_solution.bound) / ((1e-10) + abs(new_sp_solution.objective))

        new_sp_solution.add_column = true
        new_sp_solution.column = sp_core_solution[source_node].column
        new_sp_solution.column_cost = new_sp_solution.objective 
    end
    #---------------------------------------------------------------------------

    return new_sp_solution
end