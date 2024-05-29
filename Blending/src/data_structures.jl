mutable struct ModelDataNodalStochastic
    I::Int64 # number of potential input tanks
    J::Int64 # number of blending/product tanks
    T::Int64 # number of time periods

    prob::Vector{Float64} # total probability of each node in the tree

    # α::Vector{Float64} # discount factor

    # Slope/intercept of price line for customer j in scenario s and time period t is given by price_slope/intercept[s,t][j]
    price_slope::Matrix{Float64}
    price_intercept::Matrix{Float64}

    # C[i]: Capacity of input tank i
    C::Vector{Float64}

    # q[i,t]:: cost of installation for input tank i and time period t
    q::Matrix{Float64}

    # r[i,j]: physical distance/transportation cost between input tank i and blending tank j 
    r::Matrix{Float64}

    # b[i,t]: production cost at input tank i in time period t
    b::Matrix{Float64}

    # λ[i]: Product quality spec at input tank i 
    λ::Vector{Float64}

    # d_max[n][j]: Maximum demand for blending tank j at node n
    d_max::Vector{Vector{Float64}}

    # d_min[n][j]: Minimum demand for blending tank j at node n
    d_min::Vector{Vector{Float64}}
end

@kwdef mutable struct SolutionInfo
    status::Any
    objective::Any
    bound::Any
    solvetime::Any
    gap::Any
    decisions::Any
    num_variables::Union{Int64, Nothing} = nothing
    num_constraints::Union{Int64, Nothing} = nothing
    num_structural_constraints::Union{Int64, Nothing} = nothing
end

@kwdef mutable struct RMPSolutionInfo
    status::Any = nothing
    objective::Any = nothing
    bound::Any = nothing
    gap::Any = nothing
    solvetime::Any = nothing
    x_val = nothing
    rho_val = nothing
    gamma1 = nothing
    gamma2 = nothing
    mu = nothing
end

@kwdef mutable struct SPSolutionInfo
    status = nothing
    objective = nothing
    bound = nothing
    gap = nothing
    solvetime = nothing
    column::Union{Float64, Vector{Float64}, Nothing} = nothing
    column_cost::Union{Float64, Vector{Float64}, Nothing} = nothing
    add_column::Bool = false
end

@kwdef mutable struct CGSolutionInfo
    UB = nothing
    LB = nothing
    relaxed_RMP_solution::Union{RMPSolutionInfo, Nothing} = nothing
    MIP_RMP_solution::Union{RMPSolutionInfo, Nothing} = nothing
    solvetime = nothing
    solvetime_perfect_parallelization = nothing
    solution_integral::Bool = true
    iterations = nothing
    total_columns_generated = nothing
    column_generated = nothing
    columns_cost = nothing
    relative_gap_each_iter = nothing
    final_relative_gap = nothing
    initial_col_gen_time = nothing
    SP_solvetime_per_iteration = nothing
    CS_solvetime_per_iteration = nothing
    RMP_solvetime_per_iteration = nothing
    CS_iterations = nothing
end
