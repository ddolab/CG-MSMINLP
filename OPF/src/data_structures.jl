mutable struct TestNetwork
    buses::Union{UnitRange{Int64}, Vector{Int64}} # indices of buses in the network (starts from 1)
    lines::Union{UnitRange{Int64}, Vector{Int64}} # indices of lines in the network
    parent
    children
end

mutable struct TestSystemData
    buses::Union{UnitRange{Int64}, Vector{Int64}} # indices of buses in the network
    lines::Union{UnitRange{Int64}, Vector{Int64}} # indices of lines in the network
    r
    x
    A
    pd_base
    qd_base
    b
    g
    Cp_base
    C_VoLL
    pg_base
    qg_base
    v_lb
    v_ub
end

mutable struct ModelDataNodalStochastic
    B::Union{UnitRange{Int64}, Vector{Int64}} # indices of buses in the network
    L::Union{UnitRange{Int64}, Vector{Int64}} # indices of lines in the network
    parent
    children
    
    M::Int64 # number of potential gensets
    T::Int64 # number of time periods
    H::Int64 # number of hours per time period (length of scheduling horizon)

    alpha::Array{Number}

    Ct::Matrix{Float64}

    Cp::Array{Float64,3}

    C_VoLL::Array{Float64}

    #-------------- uncertain parameters --------------#
    # pd[n][i,h]
    pd::Vector{Matrix{Float64}}
    # qd[n][i,h]
    qd::Vector{Matrix{Float64}}
    #--------------------------------------------------#

    Ca::Array{Float64}

    pa::Array{Float64}

    r

    g::Array{Float64}

    x

    A

    b::Array{Float64}

    pg_ub::Array{Float64}

    qg_ub::Array{Float64}

    v_lb::Array{Float64}

    v_ub::Array{Float64}

    prob::Vector{Float64} # Total probability of a node
end

@kwdef mutable struct SolutionInfo
    status::Any = nothing
    objective::Any = nothing
    bound::Any = nothing
    solvetime::Any = nothing
    gap::Any = nothing
    decisions::Any = nothing
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
    y_val = nothing
    z_val = nothing
    rho_val = nothing
    gamma1 = nothing
    gamma2 = nothing
    gamma3 = nothing
    mu = nothing
end

@kwdef mutable struct SPColumnInfo
    add_column::Bool = false
    column = nothing
    cost::Union{Float64, Vector{Float64}, Nothing} = nothing
end

@kwdef mutable struct SPSolutionInfo
    status = nothing
    solvetime = nothing
    bound = nothing
    columns_info = nothing
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

@kwdef mutable struct ExpectedValueSolution
    status = nothing
    num_of_columns = nothing
    column = nothing
    objective = nothing
    solvetime = nothing
end
