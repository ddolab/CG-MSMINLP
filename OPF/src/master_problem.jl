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

    #-----------------Construct model---------------------
    if relaxed
        @variables(model, begin
            rho[n in 2:N, k in 1:NC[n]] >= 0
            0 <= y[m in 1:data.M, i in data.B, n in 1:(minimum(leaves(tree)[1])-1)] <= 1
            0 <= z[m in 1:data.M, n in 2:(minimum(leaves(tree)[1])-1)] <= 1
        end)
    else
        @variables(model, begin
            rho[n in 2:N, k in 1:NC[n]], Bin
            0 <= y[m in 1:data.M, i in data.B, n in 1:(minimum(leaves(tree)[1])-1)] <= 1 # forced to be binary by rho
            z[m in 1:data.M, n in 2:(minimum(leaves(tree)[1])-1)], Bin
        end)
    end

    @constraints(model, begin
        con1[m in 1:data.M, i in data.B, n in 2:N], sum(rho[n,k]*w_star[n][k][m,i] for k in 1:NC[n]) == y[m,i,tree.parent[n]]

        con2[m in 1:data.M, i in data.B, n in 2:(minimum(leaves(tree)[1])-1)], y[m,i,n]-sum(rho[n,k]*w_star[n][k][m,i] for k in 1:NC[n]) - z[m,n] <= 0
        con3[m in 1:data.M, i in data.B, n in 2:(minimum(leaves(tree)[1])-1)], z[m,n] + y[m,i,n] + sum(rho[n,k]*w_star[n][k][m,i] for k in 1:NC[n]) <= 2

        convexity[n in 2:N], sum(rho[n,k] for k in 1:NC[n]) == 1
    end)

    @objective(model, Min, sum(sum(w_cost[n][k]*rho[n,k] for k in 1:NC[n]) for n in 2:N) + 
        sum(data.prob[n]*data.Ct[m,stage(tree,n)[1]+1]*z[m,n] for m in 1:data.M, n in 2:(minimum(leaves(tree)[1])-1)))
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
    RMP_solution.y_val = value.(y)
    RMP_solution.z_val = value.(z)
    RMP_solution.rho_val = value.(rho)
    #---------------------------------------------------------------------------

    #-----------------------------RETRIEVING DUAL VALUES------------------------
    if relaxed
        if has_duals(model)
            RMP_solution.gamma1 = dual.(con1)
            RMP_solution.gamma2 = dual.(con2)
            RMP_solution.gamma3 = dual.(con3)
            RMP_solution.mu = dual.(convexity)
        else
            error("Dual solution NOT generated!")
        end
    end
    #---------------------------------------------------------------------------

    return RMP_solution
end