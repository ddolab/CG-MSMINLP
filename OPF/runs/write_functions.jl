function create_spreadsheet_FS(sol::SolutionInfo, model_type, instance_number)
    
    filename = string(model_type,"_instance_",instance_number,".xlsx")
    XLSX.openxlsx(filename, mode = "w") do xf
        sheet = xf[1]
        # Write column headers for statistics
        sheet["A1"] = "instance"
        sheet["B1"] = "obj"
        sheet["C1"] = "gap"
        sheet["D1"] = "time"
        
        # Write the statistics in row 2
        sheet["A2"] = instance_number
        sheet["B2"] = sol.objective
        sheet["C2"] = sol.gap*100 # converted to %
        sheet["D2"] = sol.solvetime
    end

end

function write_CG_gap_to_file(sol::CGSolutionInfo, model_type, instance_number)

    file = open(string(model_type, "_",instance_number, "_gap.txt"), "a")

    for gap in sol.relative_gap_each_iter
        println(file, gap*100) # converted to %
    end

    close(file)
end

function create_spreadsheet_CG(sol::CGSolutionInfo, model_type, instance_number)

    filename = string(model_type,"_instance_",instance_number,".xlsx")
    
    XLSX.openxlsx(filename, mode = "w") do xf
        sheet = xf[1]

        # Write column headers for statistics
        sheet["A1"] = "instance"
        sheet["B1"] = "obj"
        sheet["C1"] = "gapCG"
        sheet["D1"] = "gapMIP"
        sheet["E1"] = "finalGap"
        sheet["F1"] = "solvetime(s)"
        sheet["G1"] = "initial_cg_time"
        sheet["H1"] = "sp_time"
        sheet["I1"] = "cs_time"
        sheet["J1"] = "rmp_time"
        sheet["K1"] = "rmp_mip_time"
        sheet["L1"] = "cols_generated"
        sheet["M1"] = "total_iters"
        sheet["N1"] = "cs_iters"
        sheet["O1"] = "integral"
        sheet["P1"] = "solvetime_perfectly_parallel(s)"   
        
        # Write the instance number
        sheet["A2"] = instance_number

        # Write the statistics
        if sol.solution_integral
            sheet["B2"] = sol.UB
            sheet["D2"] = sol.final_relative_gap*100    # converted to %
            sheet["E2"] = sol.final_relative_gap*100    # converted to %
            sheet["K2"] = "NA"
            sheet["O2"] = "yes"
        else
            sheet["B2"] = sol.MIP_RMP_solution.objective
            sheet["D2"] = sol.MIP_RMP_solution.gap*100  # converted to %
            sheet["E2"] = 100*abs(sol.MIP_RMP_solution.objective - sol.LB) / ((1e-10) + abs(sol.MIP_RMP_solution.objective))
            sheet["K2"] = sol.MIP_RMP_solution.solvetime
            sheet["O2"] = "no"
        end

        sheet["C2"] = sol.final_relative_gap*100   # converted to %
        sheet["F2"] = sol.solvetime
        sheet["G2"] = sol.initial_col_gen_time
        sheet["H2"] = sum(sol.SP_solvetime_per_iteration)

        if length(sol.CS_solvetime_per_iteration) != 0
            sheet["I2"] = sum(sol.CS_solvetime_per_iteration)
        else
            sheet["I2"] = 0.0
        end

        sheet["J2"] = sum(sol.RMP_solvetime_per_iteration)
        sheet["L2"] = sol.total_columns_generated
        sheet["M2"] = sol.iterations
        sheet["N2"] = sol.CS_iterations
        sheet["P2"] = sol.solvetime_perfect_parallelization
    end

end