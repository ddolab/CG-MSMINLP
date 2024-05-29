function solve_stoch_model_FS(data::ModelDataNodalStochastic;
    tree::Tree,
    timelimit,
    gap
)

    N = maximum(nodes(tree))

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "NonConvex", 2)
    set_optimizer_attribute(model, "TimeLimit", timelimit)
    set_optimizer_attribute(model, "MIPGap", gap)

    # define variables
    @variables(
        model,
        begin
            x[i in 1:data.I, n in 1:(minimum(leaves(tree)[1])-1)], Bin
            data.d_min[n][j] <= d[j in 1:data.J, n in 2:N] <= data.d_max[n][j]
            0 <= c[j in 1:data.J, n in 2:N] <= 1
            F[i in 1:data.I, j in 1:data.J, n in 2:N] >= 0
        end
    )

    # capacity expansion
    @constraint(
        model,
        [i in 1:data.I, n in 2:(minimum(leaves(tree)[1])-1)], sum(x[i, n_] for n_ in root(tree, n)) <= 1
    )

    # Feasible production
    @constraint(
        model,
        [i in 1:data.I, n in 2:N],
        sum(F[i, j, n] for j in 1:data.J) <= data.C[i]*sum(x[i, n_] for n_ in root(tree, tree.parent[n]))
    )

    # demand satisfcation
    @constraint(
        model,
        [j in 1:data.J, n in 2:N],
        sum(F[i, j, n] for i in 1:data.I) == d[j, n]
    )

    # blending
    @constraint(
        model,
        [j in 1:data.J, n in 2:N],
        sum(data.λ[i] * F[i, j, n] for i in 1:data.I) == c[j, n] * d[j, n]
    )
    
    @expression(model, expr_x, sum(data.prob[n]*data.q[i,stage(tree, n)[1]+1]*x[i,n] for i in 1:data.I, n in 1:(minimum(leaves(tree)[1])-1)))
    @expression(model, expr_b, sum(data.prob[n]*data.b[i,stage(tree, n)[1]]*F[i,j,n] for j in 1:data.J, i in 1:data.I, n in 2:N ))
    @expression(model, expr_r, sum(data.prob[n]*data.r[i,j]*F[i,j,n] for i in 1:data.I, j in 1:data.J, n in 2:N))
    @expression(model, expr_F, sum(data.prob[n]*(data.price_slope[j,stage(tree, n)[1]] * c[j, n] + data.price_intercept[j,stage(tree, n)[1]]) * F[i,j,n] for i in 1:data.I, j in 1:data.J, n in 2:N))

    @objective(model, Max, expr_F - expr_x - expr_b  - expr_r)
    
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
        sol = (value.(x), value.(d), value.(c), value.(F))
    end

    num_of_vars = num_variables(model)
    num_of_constraints = num_constraints(model; count_variable_in_set_constraints = true)
    num_of_structural_constraints = num_constraints(model; count_variable_in_set_constraints = false)
    #---------------------------------------------------------------------------

    return SolutionInfo(status, objective, bound, solveTime, gap, sol, num_of_vars, num_of_constraints, num_of_structural_constraints)
end
