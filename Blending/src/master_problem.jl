# Args: data - Model data
#       w_cost - column costs
#       w_star - columns
#       NC - number of columns (a vector)
#       tree - Scenario tree
#       relaxed - whether the model is relaxed
# Returns: RMP_solution
function solve_RMP(data::ModelDataNodalStochastic, 
    w_cost, 
    w_star, 
    NC;
    tree, 
    relaxed = true)

    N = maximum(nodes(tree))

    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    set_optimizer_attribute(model, "TimeLimit", 7200.0)

    if relaxed
        set_optimizer_attribute(model, "Presolve", 0)
        set_optimizer_attribute(model, "Method", 2)
    end
    #--------------------------------------------

    #-----------------Construct model---------------------
    if relaxed
        @variables(model, begin
            rho[n in 2:N, k in 1:NC[n]] >= 0
            0 <= x[i in 1:data.I, n in 1:(minimum(leaves(tree)[1])-1)] <= 1
        end)
    else
        @variables(model, begin
            rho[n in 2:N, k in 1:NC[n]], Bin
            x[i in 1:data.I, n in 1:(minimum(leaves(tree)[1])-1)], Bin
        end)
    end

    @constraints(model, begin
        con1[i in 1:data.I, n in 2:(minimum(leaves(tree)[1])-1)], sum(rho[n,k]*w_star[n][k,i] for k in 1:NC[n]) + x[i,n] <= 1
        con2[i in 1:data.I, n in 2:N], sum(rho[n,k]*w_star[n][k,i] for k in 1:NC[n]) == sum(x[i, n_] for n_ in root(tree, tree.parent[n]))
        convexity[n in 2:N], sum(rho[n,k] for k in 1:NC[n]) == 1
    end)

    @objective(model, Min, sum(data.prob[n]*data.q[i,stage(tree, n)[1]+1]*x[i,n] for i in 1:data.I, n in 1:(minimum(leaves(tree)[1])-1) ) +
        sum(sum(w_cost[n][k]*rho[n,k] for k in 1:NC[n]) for n in 2:N) )
    #---------------------------------------------------------------------------

    optimize!(model)


    #-----------------------------RETRIEVING SOLUTION---------------------------
    RMP_solution = RMPSolutionInfo()
    RMP_solution.status = termination_status(model)
    RMP_solution.solvetime = solve_time(model)

    RMP_solution.objective = objective_value(model)
    if !relaxed
        RMP_solution.bound = objective_bound(model)
        RMP_solution.gap = abs(RMP_solution.objective - RMP_solution.bound) / ((1e-10) + abs(RMP_solution.objective))
    end
    #---------------------------------------------------------------------------

    #--------------------------Storing variable values--------------------------
    RMP_solution.x_val = value.(x)
    RMP_solution.rho_val = value.(rho)
    #---------------------------------------------------------------------------

    #-----------------------------RETRIEVING DUAL VALUES------------------------
    if relaxed
        if has_duals(model)
            RMP_solution.gamma1 = dual.(con1)
            RMP_solution.gamma2 = dual.(con2)
            RMP_solution.mu = dual.(convexity)
        else
            error("Dual solution NOT generated!")
        end
    end
    #---------------------------------------------------------------------------

    return RMP_solution
end