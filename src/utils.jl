function value_to_df_2dim(var)
    solution = DataFrame(var.data, :auto)
    ax1 = var.axes[1]
    ax2 = var.axes[2]
    cols = names(solution)
    insertcols!(solution, 1, :r_id => ax1)
    solution = stack(solution, Not(:r_id), variable_name=:hour)
    solution.hour = foldl(replace, [cols[i] => ax2[i] for i in 1:length(ax2)], init=solution.hour)
    rename!(solution, :value => :gen)
    solution.hour = convert.(Int64,solution.hour)
    return solution
end

function process_data(path)
    datadir = joinpath(path) 

    gen_info = CSV.read(joinpath(datadir,"Generators_data.csv"), DataFrame);
    fuels = CSV.read(joinpath(datadir,"Fuels_data.csv"), DataFrame);
    loads = CSV.read(joinpath(datadir,"Load_data.csv"), DataFrame);
    gen_var = CSV.read(joinpath(datadir,"Generators_variability.csv"), DataFrame);
    network = CSV.read(joinpath(datadir,"Network.csv"), DataFrame);

    # Rename all columns to lowercase (by convention)
    for f in [gen_info, fuels, loads, gen_var, network]
        rename!(f,lowercase.(names(f)))
    end

    # remove irrelevant gen_info columns
    select!(gen_info, Not(["commit", "heat_rate_mmbtu_mwh_iqr", "heat_rate_mmbtu_mwh_std"]))

    # process loads data frame
    select!(loads, Not(1:4)) # remove unnecessary columns
    rename!(loads, :time_index => :hour)
    # rename load columns to just zone numbers
    for name in names(loads)
        if occursin("load_mw_z", name)
            zone_num = replace(name, "load_mw_z" => "")
            rename!(loads, name => zone_num)
        end
    end

    # convert loads to long
    loads_long = stack(loads, Not(:hour), variable_name=:zone, value_name=:demand)
    loads_long.zone = parse.(Int, loads_long.zone)

    # create gen_df
    gen_df = outerjoin(gen_info,  fuels, on = :fuel) # load in fuel costs and add to data frame
    rename!(gen_df, :cost_per_mmbtu => :fuel_cost)   # rename column for fuel cost
    gen_df.fuel_cost[ismissing.(gen_df[:,:fuel_cost])] .= 0
    dropmissing!(gen_df)
    gen_df.fuel[gen_df.hydro .== 1] .= "hydro"

    # make gen_var columns r_id
    rename!(gen_var, names(gen_var)[1] => :hour)
    for i in 1:length(names(gen_var))-1
        rename!(gen_var, names(gen_var)[i+1] => Symbol(i))
    end
    # convert gen_var to long
    gen_var_long = stack(gen_var, 
                        Not(:hour), variable_name=:r_id,
                                            value_name=:cf)
    gen_var_long.r_id = parse.(Int, gen_var_long.r_id); # r_ids should be int

    for name in names(network)
        if occursin("z", name)
            zone_num = replace(name, "z" => "")
            rename!(network, name => zone_num)
        end
    end

    return gen_df, gen_var_long, loads_long, network
end

function evaluate_commitment_on_actuals(
    uc_solution,
    gen_df,
    loads_actual,
    gen_var_actual,
    network,
    mip_gap;
    cost_NSE = 3000.0
)
    """
    Fix commitment decisions from day-ahead UC, 
    re-dispatch generation/storage/reserves against actual conditions
    """
    
    UC = Model(HiGHS.Optimizer)
    set_optimizer_attribute(UC, "mip_rel_gap", mip_gap)
    
    # Same setup as UC model
    storage_df = gen_df[gen_df.stor .== 1, :]
    gen_df = gen_df[gen_df.stor .== 0, :]
    
    G_thermal = gen_df[gen_df[!,:up_time] .> 0,:r_id]
    G_nonthermal = gen_df[gen_df[!,:up_time] .== 0,:r_id]       
    G_var = gen_df[gen_df[!,:vre] .== 1, :r_id]
    G_nonvar = gen_df[gen_df[!,:vre] .== 0,:r_id]
    G_nt_nonvar = intersect(G_nonvar, G_nonthermal)
    B = storage_df.r_id
    G = gen_df.r_id
    L = network.network_lines
    T = unique(loads_actual.hour)
    T_red = T[1:end-1]
    Z = unique(gen_df.zone)
    
    # Define all variables
    @variables(UC, begin
        GEN[i in G, t in T] >= 0
        GENAUX[i in G_thermal, t in T] >= 0
        CHARGE[b in B, t in T] >= 0
        DISCHARGE[b in B, t in T] >= 0
        SOC[b in B, t in T] >= 0
        FLOW[l in L, t in T]
        RESUP[i in G_thermal, t in T] >= 0
        RESDN[i in G_thermal, t in T] >= 0
        # Load shedding as slack variable
        LOADSHED[z in Z, t in T] >= 0
    end)
    
    # **KEY: Fix commitment from day-ahead solution**
    @variable(UC, COMMIT[i in G_thermal, t in T], Bin)
    
    # Extract commitment values - handle column naming
    commit_df = uc_solution.commit
    for i in G_thermal
        for t in T
            commit_val = commit_df[(commit_df.r_id .== i) .& (commit_df.hour .== t), :gen]
            # if length(commit_val) > 0
            fix(COMMIT[i, t], commit_val[1]; force=true)
            # end
        end
    end
    # Objective: minimize redispatch cost + load shedding penalty
    @objective(UC, Min,
        sum(gen_df[gen_df.r_id .== i, :var_om_cost_per_mwh][1] * GEN[i, t] 
            for i in G, t in T) +
        sum(cost_NSE * LOADSHED[z, t] for z in Z, t in T)
    )
    
    # Network balance with load shedding
    @constraint(UC, NodalBalance[z in Z, t in T],
        sum(GEN[i, t] for i in G if gen_df[gen_df.r_id .== i, :zone][1] == z) +
        sum(DISCHARGE[b, t] - CHARGE[b, t] for b in B if storage_df[storage_df.r_id .== b, :zone][1] == z) +
        LOADSHED[z, t] -
        loads_actual[loads_actual.zone .== z .&& loads_actual.hour .== t, :demand][1] -
        sum(network[l, Symbol("$z")] * FLOW[l, t] for l in L) == 0
    )
    
    # Line flow limits
    @constraint(UC, FlowMax[l in L, t in T],
        FLOW[l, t] <= network[l, :line_max_flow_mw]
    )
    @constraint(UC, FlowMin[l in L, t in T],
        FLOW[l, t] >= -network[l, :line_max_flow_mw]
    )
    
    # Thermal capacity limits (respecting fixed commitment)
    @constraint(UC, Cap_thermal_min[i in G_thermal, t in T],
        GEN[i, t] >= gen_df[gen_df.r_id .== i, :min_power][1] * 
                     gen_df[gen_df.r_id .== i, :existing_cap_mw][1] * COMMIT[i, t]
    )
    @constraint(UC, Cap_thermal_max[i in G_thermal, t in T],
        GEN[i, t] <= gen_df[gen_df.r_id .== i, :existing_cap_mw][1] * COMMIT[i, t]
    )
    
    # Non-thermal non-variable capacity
    @constraint(UC, Cap_nt_nonvar[i in G_nt_nonvar, t in T],
        GEN[i, t] <= gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    
    # Variable generation using ACTUAL capacity factors
    @constraint(UC, Cap_var[i in G_var, t in T],
        GEN[i, t] <= gen_var_actual[(gen_var_actual.r_id .== i) .& (gen_var_actual.hour .== t), :cf][1] * 
                     gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    
    # Battery constraints
    @constraint(UC, ChargeCap[b in B, t in T],
        CHARGE[b, t] <= storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
    )
    @constraint(UC, DischargeCap[b in B, t in T],
        DISCHARGE[b, t] <= storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
    )
    @constraint(UC, SOCCap[b in B, t in T],
        SOC[b, t] <= storage_df[storage_df.r_id .== b, :existing_cap_mwh][1]
    )
    
    @constraint(UC, cStateOfChargeEnd[b in B],
        SOC[b, T[end]] == storage_df[storage_df.r_id .== b, :existing_cap_mwh][1] * 0.5
    )
    @constraint(UC, cStateOfChargeStart[b in B],
        SOC[b, T[1]] == storage_df[storage_df.r_id .== b, :existing_cap_mwh][1] * 0.5 +
        storage_df[storage_df.r_id .== b, :charge_efficiency][1] * CHARGE[b, T[1]] -
        DISCHARGE[b, T[1]] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1]
    )
    @constraint(UC, cStateOfCharge[b in B, t in T[2:end]],
        SOC[b, t] == SOC[b, t-1] +
        storage_df[storage_df.r_id .== b, :charge_efficiency][1] * CHARGE[b, t] -
        DISCHARGE[b, t] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1]
    )
    
    # Auxiliary generation variable for ramping
    @constraint(UC, AuxGen[i in G_thermal, t in T],
        GENAUX[i, t] == GEN[i, t] - 
        gen_df[gen_df.r_id .== i, :min_power][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1] * COMMIT[i, t]
    )
    
    # Ramping constraints (thermal)
    @constraint(UC, RampUp_thermal[i in G_thermal, t in T_red],
        GENAUX[i, t+1] - GENAUX[i, t] <= 
        gen_df[gen_df.r_id .== i, :ramp_up_percentage][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    @constraint(UC, RampDn_thermal[i in G_thermal, t in T_red],
        GENAUX[i, t] - GENAUX[i, t+1] <= 
        gen_df[gen_df.r_id .== i, :ramp_dn_percentage][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    
    # Ramping constraints (non-thermal)
    @constraint(UC, RampUp_nonthermal[i in G_nonthermal, t in T_red],
        GEN[i, t+1] - GEN[i, t] <= 
        gen_df[gen_df.r_id .== i, :ramp_up_percentage][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    @constraint(UC, RampDn_nonthermal[i in G_nonthermal, t in T_red],
        GEN[i, t] - GEN[i, t+1] <= 
        gen_df[gen_df.r_id .== i, :ramp_dn_percentage][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    
    # Reserve constraints
    @constraint(UC, ResUpCap[i in G_thermal, t in T],
        RESUP[i, t] <= gen_df[gen_df.r_id .== i, :existing_cap_mw][1] * COMMIT[i, t] - GEN[i, t]
    )
    @constraint(UC, ResDnCap[i in G_thermal, t in T],
        RESDN[i, t] <= GEN[i, t] - 
        gen_df[gen_df.r_id .== i, :min_power][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1] * COMMIT[i, t]
    )
    
    @constraint(UC, ResUpRamp[i in G_thermal, t in T],
        RESUP[i, t] <= gen_df[gen_df.r_id .== i, :ramp_up_percentage][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    @constraint(UC, ResDnRamp[i in G_thermal, t in T],
        RESDN[i, t] <= gen_df[gen_df.r_id .== i, :ramp_dn_percentage][1] * 
        gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
    )
    
    ResReqUp = Dict(t => 300 for t in T)
    ResReqDn = Dict(t => 0 for t in T)
    
    @constraint(UC, ResUpRequirement[t in T],
        sum(RESUP[i, t] for i in G_thermal) >= ResReqUp[t]
    )
    @constraint(UC, ResDnRequirement[t in T],
        sum(RESDN[i, t] for i in G_thermal) >= ResReqDn[t]
    )
    
    optimize!(UC)
    
    # Extract results
    gen_actual = value_to_df_2dim(value.(GEN))
    loadshed = value_to_df_2dim(value.(LOADSHED))
    resup = value_to_df_2dim(value.(RESUP))
    resdn = value_to_df_2dim(value.(RESDN))
    
    # Calculate curtailment
    curtail = DataFrame(r_id = Int[], hour = Int[], curt = Float64[])
    for i in G_var
        for t in T
            available = gen_var_actual[(gen_var_actual.r_id .== i) .& (gen_var_actual.hour .== t), :cf][1] * 
                       gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
            actual_gen = value(GEN[i, t])
            curt_val = max(0, available - actual_gen)
            push!(curtail, (i, t, curt_val))
        end
    end
    
    # Calculate total CO2 emissions
    total_co2 = sum(gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    value(GEN[i,t]) for i in G_nonvar for t in T) +
                sum(gen_df[gen_df.r_id .== i,:start_fuel_mmbtu_per_mw][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                    value(COMMIT[i,T[1]]) for i in G_thermal)

    # battery vars
    charge_vals = Dict(b => [value(CHARGE[b, t]) for t in T] for b in B)
    discharge_vals = Dict(b => [value(DISCHARGE[b, t]) for t in T] for b in B)

    charge_df = DataFrame(
        r_id = repeat(B, inner=length(T)),
        hour = repeat(collect(T), outer=length(B)),
        value = vcat([charge_vals[b] for b in B]...)
    )

    discharge_df = DataFrame(
        r_id = repeat(B, inner=length(T)),
        hour = repeat(collect(T), outer=length(B)),
        value = vcat([discharge_vals[b] for b in B]...)
    )
    # Extract flow values into a DataFrame
    flow_vals = Dict((l, t) => value(FLOW[l, t]) for l in L, t in T)

    flows = DataFrame(
        line = repeat(L, inner=length(T)),
        hour = repeat(collect(T), outer=length(L)),
        flow = [flow_vals[(l, t)] for l in L for t in T]
    )
    
    return (
        model = UC,
        gen = gen_actual,
        commit = uc_solution.commit,  # Keep original commitments
        loadshed = loadshed,
        curtail = curtail,
        reserves_up = resup,
        reserves_down = resdn,
        charge_df = charge_df,
        discharge_df = discharge_df,
        flows = flows,
        total_loadshed = sum(value.(LOADSHED)),
        total_curtailment = sum(curtail.curt),
        total_co2 = total_co2,
        objective = objective_value(UC)
    )
end