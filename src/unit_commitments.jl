function unit_commitment_no_transport(gen_df, loads, gen_var, mip_gap)
    UC = Model(HiGHS.Optimizer)
    set_optimizer_attribute(UC, "mip_rel_gap", mip_gap) 

    # isolate storage
    storage_df = gen_df[gen_df.stor .== 1, :]
    gen_df = gen_df[gen_df.stor .== 0, :]

    # gen subsets
    G_thermal = gen_df[gen_df[!,:up_time] .> 0,:r_id] # does NOT correspond to therm in gen_df
    G_nonthermal = gen_df[gen_df[!,:up_time] .== 0,:r_id]       
    G_var = gen_df[gen_df[!,:vre] .== 1, :r_id]
    G_nonvar = gen_df[gen_df[!,:vre] .== 0,:r_id]
    G_nt_nonvar = intersect(G_nonvar, G_nonthermal)
    B = storage_df.r_id # batteries
    G = gen_df.r_id
    T = loads.hour
    T_red = loads.hour[1:end-1]  # reduced time periods without last one

    # Decision variables   
    @variables(UC, begin
        GEN[G, T]  >= 0     # generation
        GENAUX[G_thermal, T] >= 0 # auxiliary generation variable
        COMMIT[G_thermal, T], Bin # commitment status (Bin=binary)
        START[G_thermal, T], Bin  # startup decision
        SHUT[G_thermal, T], Bin   # shutdown decision
        RESUP[G_thermal, T]       # up reserve capacity
        RESDN[G_thermal, T]       # down reserve capacity

        # battery vars
        0 <= SOC[b in B, T] <= storage_df[storage_df.r_id .== b, :existing_cap_mw][1]  # state of charge
        0 <= CHARGE[b in B, T] <= storage_df[storage_df.r_id .== b, :ramp_up_percentage][1] * 
                                    storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
        0 <= DISCHARGE[b in B, T] <= storage_df[storage_df.r_id .== b, :ramp_dn_percentage][1] * 
                                    storage_df[storage_df.r_id .== b, :existing_cap_mw][1]  
    end)
                
    # Objective function
    @objective(UC, Min, 
        sum( (gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * gen_df[gen_df.r_id .== i,:fuel_cost][1] +
            gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1]) * GEN[i,t] 
                        for i in G_nonvar for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1] * GEN[i,t] 
                        for i in G_var for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:start_cost_per_mw][1] * 
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
            START[i,t] 
                        for i in G_thermal for t in T)
    )
    
    # Demand constraint
    @constraint(UC, cDemand[t in T], 
        sum(GEN[i,t] for i in G) + sum(DISCHARGE[b,t] - CHARGE[b,t] for b in B) == loads[loads.hour .== t,:demand][1])

    # Capacity constraints (thermal generators requiring commitment)
    @constraint(UC, Cap_thermal_min[i in G_thermal, t in T], 
        GEN[i,t] >= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    @constraint(UC, Cap_thermal_max[i in G_thermal, t in T], 
        GEN[i,t] <= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # Capacity constraints (non-variable generation not requiring commitment)
    @constraint(UC, Cap_nt_nonvar[i in G_nt_nonvar, t in T], 
        GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # Capacity constraints (variable generation)

    @constraint(UC, Cap_var[i in G_var, t in T], 
        GEN[i,t] <= gen_var[(gen_var.hour .== t) .& (gen_var.r_id .== i), :cf][1] *
                    gen_df[gen_df.r_id .== i, :existing_cap_mw][1])
    
    # battery constraints 
    @constraint(UC, cStateOfChargeEnd[b in B],
        SOC[b, T[end]] == 0.5 * storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
    )
    @constraint(UC, cStateOfChargeStart[b in B],
        SOC[b, T[1]] == 0.5 * storage_df[storage_df.r_id .== b, :existing_cap_mw][1] + # start charge
                (CHARGE[b, T[1]] * storage_df[storage_df.r_id .== b, :charge_efficiency][1] -
                DISCHARGE[b, T[1]] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1])
    )
    @constraint(UC, cStateOfCharge[b in B, t in T[2:end]],
        SOC[b, t] == SOC[b, t-1] +
                    (CHARGE[b, t] * storage_df[storage_df.r_id .== b, :charge_efficiency][1] -
                    DISCHARGE[b, t] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1])
    )

    # Unit commitment constraints
    @constraint(UC, Startup[i in G_thermal, t in T],
        COMMIT[i, t] >= sum(START[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:up_time][1]):t)))

    @constraint(UC, Shutdown[i in G_thermal, t in T],
        1-COMMIT[i, t] >= sum(SHUT[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:down_time][1]):t)))
    
    @constraint(UC, CommitmentStatus[i in G_thermal, t in T_red],
        COMMIT[i,t+1] - COMMIT[i,t] == START[i,t+1] - SHUT[i,t+1])
    
    # New auxiliary variable GENAUX for generation above the minimum output level
    # for committed thermal units (only created for thermal generators)
    @constraint(UC, AuxGen[i in G_thermal, t in T],
        GENAUX[i,t] == GEN[i,t] - COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    
    # Ramp equations for thermal generators (constraining GENAUX)
    @constraint(UC, RampUp_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t+1] - GENAUX[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] )

    @constraint(UC, RampDn_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t] - GENAUX[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] )

    # Ramp equations for non-thermal generators (constraining total generation GEN)
    @constraint(UC, RampUp_nonthermal[i in G_nonthermal, t in T_red], 
        GEN[i,t+1] - GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] )

    @constraint(UC, RampDn[i in G, t in T_red], 
        GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] )
    
                                 
    # Reserve equations
    # (1) Reserves limited by committed capacity of generator
    @constraint(UC, ResUpCap[i in G_thermal, t in T],
        RESUP[i,t] <= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1]
                        - GEN[i,t])

    @constraint(UC, ResDnCap[i in G_thermal, t in T],
        RESDN[i,t] <= GEN[i,t] - 
                        COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])

    # (2) Reserves limited by ramp rates
    @constraint(UC, ResUpRamp[i in G_thermal, t in T],
        RESUP[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                          gen_df[gen_df.r_id .== i,:ramp_up_percentage][1])

    @constraint(UC, ResDnRamp[i in G_thermal, t in T],
        RESDN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                          gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])

    # Compute reserve requirements
    ResReqUp = Dict(T[i] => 300  
              + 0.05 * loads[loads.hour .== T[i],:demand][1]
                    for i in 1:length(T))
    ResReqDn = Dict(T[i] => 0
              + 0.05 * loads[loads.hour .== T[i],:demand][1]                    
                    for i in 1:length(T))
    # (3) Overall reserve requirements
    @constraint(UC, ResUpRequirement[t in T],
        sum(RESUP[i,t] for i in G_thermal) >= 
                ResReqUp[t])

    @constraint(UC, ResDnRequirement[t in T],
        sum(RESDN[i,t] for i in G_thermal) >= 
                ResReqDn[t])
    
    # Solve statement (! indicates runs in place)
    optimize!(UC)

    # Generation solution
    gen = value_to_df_2dim(value.(GEN))

    # Commitment status solution
    commit = value_to_df_2dim(value.(COMMIT))

    # Calculate curtailment
    curtail = DataFrame(r_id = Int[], hour = Int[], curt = Float64[])
    for i in G_var
        for t in T
            cf = gen_var[(gen_var.hour .== t) .& (gen_var.r_id .== i), :cf][1]
            existing_cap = gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
            generation = value(GEN[i,t])
            curt_val = cf * existing_cap - generation
            push!(curtail, (r_id = i, hour = t, curt = curt_val))
        end
    end

    # Calculate total CO2 emissions
    total_co2 = sum(gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    value(GEN[i,t]) for i in G_nonvar for t in T) +
                sum(gen_df[gen_df.r_id .== i,:start_fuel_mmbtu_per_mw][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                    value(START[i,t]) for i in G_thermal for t in T)

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

    # Return the solution and objective as named tuple
    return (
        gen,
        commit,
        curtail,
        total_co2,
        charge_df,
        discharge_df,
        reserves = value_to_df_2dim(value.(RESUP)),
        cost = objective_value(UC),
        status = termination_status(UC)
    )

end

function unit_commitment_transport(gen_df, loads, gen_var, network, mip_gap)
    print(typeof(network[1,4]))
    UC = Model(HiGHS.Optimizer)
    set_optimizer_attribute(UC, "mip_rel_gap", mip_gap) 

    # isolate storage
    storage_df = gen_df[gen_df.stor .== 1, :]
    gen_df = gen_df[gen_df.stor .== 0, :]

    # gen subsets
    G_thermal = gen_df[gen_df[!,:up_time] .> 0,:r_id] # does NOT correspond to therm in gen_df
    G_nonthermal = gen_df[gen_df[!,:up_time] .== 0,:r_id]       
    G_var = gen_df[gen_df[!,:vre] .== 1, :r_id]
    G_nonvar = gen_df[gen_df[!,:vre] .== 0,:r_id]
    G_nt_nonvar = intersect(G_nonvar, G_nonthermal)
    B = storage_df.r_id # batteries
    G = gen_df.r_id
    L = network.network_lines
    T = unique(loads.hour)
    T_red = T[1:end-1]  # reduced time periods without last one
    Z = unique(gen_df.zone)

    # Decision variables   
    @variables(UC, begin
        GEN[G, T]  >= 0     # generation
        GENAUX[G_thermal, T] >= 0 # auxiliary generation variable
        COMMIT[G_thermal, T], Bin # commitment status (Bin=binary)
        START[G_thermal, T], Bin  # startup decision
        SHUT[G_thermal, T], Bin   # shutdown decision
        RESUP[G_thermal, T]       # up reserve capacity
        RESDN[G_thermal, T]       # down reserve capacity
        FLOW[l in L, t in T] >= 0 # power flows

        # battery vars
        0 <= SOC[b in B, T] <= storage_df[storage_df.r_id .== b, :existing_cap_mw][1]  # state of charge
        0 <= CHARGE[b in B, T] <= storage_df[storage_df.r_id .== b, :ramp_up_percentage][1] * 
                                    storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
        0 <= DISCHARGE[b in B, T] <= storage_df[storage_df.r_id .== b, :ramp_dn_percentage][1] * 
                                    storage_df[storage_df.r_id .== b, :existing_cap_mw][1]  
    end)
                
    # Objective function
    @objective(UC, Min, 
        sum( (gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * gen_df[gen_df.r_id .== i,:fuel_cost][1] +
            gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1]) * GEN[i,t] 
                        for i in G_nonvar for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1] * GEN[i,t] 
                        for i in G_var for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:start_cost_per_mw][1] * 
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
            START[i,t] 
                        for i in G_thermal for t in T)
    )
    
    # OPF constraints

    # Auxiliary variable: total generation per zone and time
    @variable(UC, ZONE_GEN[z in Z, t in T] >= 0)

    @constraint(UC, ZoneGenDef[z in Z, t in T],
        ZONE_GEN[z, t] == sum(GEN[i, t] for i in G if gen_df[gen_df.r_id .== i, :zone][1] == z)
    )

    # power line capacities
    @constraint(UC, FlowMax[l in L, t in T],
        FLOW[l, t] <= network.line_max_flow_mw[l]
    )

    @constraint(UC, FlowMin[l in L, t in T],
        FLOW[l, t] >= -network.line_max_flow_mw[l]
    )

    # demand balance by zone
    @constraint(UC, NodalBalance[z in Z, t in T],
        ZONE_GEN[z, t] + 
        sum(DISCHARGE[b,t] - CHARGE[b,t] for b in B if storage_df[storage_df.r_id .== b, :zone][1] == z) +
        sum(FLOW[l, t] * network[l, "$z"] * (1 - network.line_loss[l]/2) for l in L if network[l, "$z"] != 0) ==
        loads[(loads.hour .== t) .& (loads.zone .== z), :demand][1]
    )

    # Capacity constraints (thermal generators requiring commitment)
    @constraint(UC, Cap_thermal_min[i in G_thermal, t in T], 
        GEN[i,t] >= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    @constraint(UC, Cap_thermal_max[i in G_thermal, t in T], 
        GEN[i,t] <= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # Capacity constraints (non-variable generation not requiring commitment)
    @constraint(UC, Cap_nt_nonvar[i in G_nt_nonvar, t in T], 
        GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # Capacity constraints (variable generation)

    @constraint(UC, Cap_var[i in G_var, t in T], 
        GEN[i,t] <= gen_var[(gen_var.hour .== t) .& (gen_var.r_id .== i), :cf][1] *
                    gen_df[gen_df.r_id .== i, :existing_cap_mw][1])
    
    # battery constraints 
    @constraint(UC, cStateOfChargeEnd[b in B],
        SOC[b, T[end]] == 0.5 * storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
    )
    @constraint(UC, cStateOfChargeStart[b in B],
        SOC[b, T[1]] == 0.5 * storage_df[storage_df.r_id .== b, :existing_cap_mw][1] + # start charge
                (CHARGE[b, T[1]] * storage_df[storage_df.r_id .== b, :charge_efficiency][1] -
                DISCHARGE[b, T[1]] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1])
    )
    @constraint(UC, cStateOfCharge[b in B, t in T[2:end]],
        SOC[b, t] == SOC[b, t-1] +
                    (CHARGE[b, t] * storage_df[storage_df.r_id .== b, :charge_efficiency][1] -
                    DISCHARGE[b, t] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1])
    )

    # Unit commitment constraints
    @constraint(UC, Startup[i in G_thermal, t in T],
        COMMIT[i, t] >= sum(START[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:up_time][1]):t)))

    @constraint(UC, Shutdown[i in G_thermal, t in T],
        1-COMMIT[i, t] >= sum(SHUT[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:down_time][1]):t)))
    
    @constraint(UC, CommitmentStatus[i in G_thermal, t in T_red],
        COMMIT[i, t+1] - COMMIT[i, t] == START[i, t+1] - SHUT[i, t+1]
    )
    # New auxiliary variable GENAUX for generation above the minimum output level
    # for committed thermal units (only created for thermal generators)
    @constraint(UC, AuxGen[i in G_thermal, t in T],
        GENAUX[i,t] == GEN[i,t] - COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    
    # Ramp equations for thermal generators (constraining GENAUX)
    @constraint(UC, RampUp_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t+1] - GENAUX[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] )

    @constraint(UC, RampDn_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t] - GENAUX[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] )

    # Ramp equations for non-thermal generators (constraining total generation GEN)
    @constraint(UC, RampUp_nonthermal[i in G_nonthermal, t in T_red], 
        GEN[i,t+1] - GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] )

    @constraint(UC, RampDn[i in G, t in T_red], 
        GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] )
    
                                 
    # Reserve equations
    # (1) Reserves limited by committed capacity of generator
    @constraint(UC, ResUpCap[i in G_thermal, t in T],
        RESUP[i,t] <= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1]
                        - GEN[i,t])

    @constraint(UC, ResDnCap[i in G_thermal, t in T],
        RESDN[i,t] <= GEN[i,t] - 
                        COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])

    # (2) Reserves limited by ramp rates
    @constraint(UC, ResUpRamp[i in G_thermal, t in T],
        RESUP[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                          gen_df[gen_df.r_id .== i,:ramp_up_percentage][1])

    @constraint(UC, ResDnRamp[i in G_thermal, t in T],
        RESDN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                          gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])

    # Compute reserve requirements
    ResReqUp = Dict(T[i] => 300  
              + 0.05 * loads[loads.hour .== T[i],:demand][1]
                    for i in 1:length(T))
    ResReqDn = Dict(T[i] => 0
              + 0.05 * loads[loads.hour .== T[i],:demand][1]                    
                    for i in 1:length(T))
                        
    # (3) Overall reserve requirements
    @constraint(UC, ResUpRequirement[t in T],
        sum(RESUP[i,t] for i in G_thermal) >= 
                ResReqUp[t])

    @constraint(UC, ResDnRequirement[t in T],
        sum(RESDN[i,t] for i in G_thermal) >= 
                ResReqDn[t])
    
    # Solve statement (! indicates runs in place)
    optimize!(UC)

    # Generation solution
    gen = value_to_df_2dim(value.(GEN))

    # Commitment status solution
    commit = value_to_df_2dim(value.(COMMIT))

    # Calculate curtailment
    curtail = DataFrame(r_id = Int[], hour = Int[], curt = Float64[])
    for i in G_var
        for t in T
            cf = gen_var[(gen_var.hour .== t) .& (gen_var.r_id .== i), :cf][1]
            existing_cap = gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
            generation = value(GEN[i,t])
            curt_val = cf * existing_cap - generation
            push!(curtail, (r_id = i, hour = t, curt = curt_val))
        end
    end

    # Calculate total CO2 emissions
    total_co2 = sum(gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    value(GEN[i,t]) for i in G_nonvar for t in T) +
                sum(gen_df[gen_df.r_id .== i,:start_fuel_mmbtu_per_mw][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                    value(START[i,t]) for i in G_thermal for t in T)

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

    # Return the solution and objective as named tuple
    return (
        gen,
        commit,
        curtail,
        total_co2,
        charge_df,
        discharge_df,
        flows,
        reserves = value_to_df_2dim(value.(RESUP)),
        cost = objective_value(UC),
        status = termination_status(UC)
    )

end

function unit_commitment_transport(gen_df, loads, gen_var, network, mip_gap)
    print(typeof(network[1,4]))
    UC = Model(HiGHS.Optimizer)
    set_optimizer_attribute(UC, "mip_rel_gap", mip_gap) 

    # isolate storage
    storage_df = gen_df[gen_df.stor .== 1, :]
    gen_df = gen_df[gen_df.stor .== 0, :]

    # gen subsets
    G_thermal = gen_df[gen_df[!,:up_time] .> 0,:r_id] # does NOT correspond to therm in gen_df
    G_nonthermal = gen_df[gen_df[!,:up_time] .== 0,:r_id]       
    G_var = gen_df[gen_df[!,:vre] .== 1, :r_id]
    G_nonvar = gen_df[gen_df[!,:vre] .== 0,:r_id]
    G_nt_nonvar = intersect(G_nonvar, G_nonthermal)
    B = storage_df.r_id # batteries
    G = gen_df.r_id
    L = network.network_lines
    T = unique(loads.hour)
    T_red = T[1:end-1]  # reduced time periods without last one
    Z = unique(gen_df.zone)

    # Decision variables   
    @variables(UC, begin
        GEN[G, T]  >= 0     # generation
        GENAUX[G_thermal, T] >= 0 # auxiliary generation variable
        COMMIT[G_thermal, T], Bin # commitment status (Bin=binary)
        START[G_thermal, T], Bin  # startup decision
        SHUT[G_thermal, T], Bin   # shutdown decision
        RESUP[G_thermal, T]       # up reserve capacity
        RESDN[G_thermal, T]       # down reserve capacity
        FLOW[l in L, t in T] >= 0 # power flows

        # battery vars
        0 <= SOC[b in B, T] <= storage_df[storage_df.r_id .== b, :existing_cap_mw][1]  # state of charge
        0 <= CHARGE[b in B, T] <= storage_df[storage_df.r_id .== b, :ramp_up_percentage][1] * 
                                    storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
        0 <= DISCHARGE[b in B, T] <= storage_df[storage_df.r_id .== b, :ramp_dn_percentage][1] * 
                                    storage_df[storage_df.r_id .== b, :existing_cap_mw][1]  
    end)
                
    # Objective function
    @objective(UC, Min, 
        sum( (gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * gen_df[gen_df.r_id .== i,:fuel_cost][1] +
            gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1]) * GEN[i,t] 
                        for i in G_nonvar for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1] * GEN[i,t] 
                        for i in G_var for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:start_cost_per_mw][1] * 
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
            START[i,t] 
                        for i in G_thermal for t in T)
    )
    
    # OPF constraints

    # Auxiliary variable: total generation per zone and time
    @variable(UC, ZONE_GEN[z in Z, t in T] >= 0)

    @constraint(UC, ZoneGenDef[z in Z, t in T],
        ZONE_GEN[z, t] == sum(GEN[i, t] for i in G if gen_df[gen_df.r_id .== i, :zone][1] == z)
    )

    # power line capacities
    @constraint(UC, FlowMax[l in L, t in T],
        FLOW[l, t] <= network.line_max_flow_mw[l]
    )

    @constraint(UC, FlowMin[l in L, t in T],
        FLOW[l, t] >= -network.line_max_flow_mw[l]
    )

    # demand balance by zone
    @constraint(UC, NodalBalance[z in Z, t in T],
        ZONE_GEN[z, t] + 
        sum(DISCHARGE[b,t] - CHARGE[b,t] for b in B if storage_df[storage_df.r_id .== b, :zone][1] == z) +
        sum(FLOW[l, t] * network[l, "$z"] * (1 - network.line_loss[l]/2) for l in L if network[l, "$z"] != 0) ==
        loads[(loads.hour .== t) .& (loads.zone .== z), :demand][1]
    )

    # Capacity constraints (thermal generators requiring commitment)
    @constraint(UC, Cap_thermal_min[i in G_thermal, t in T], 
        GEN[i,t] >= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    @constraint(UC, Cap_thermal_max[i in G_thermal, t in T], 
        GEN[i,t] <= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # Capacity constraints (non-variable generation not requiring commitment)
    @constraint(UC, Cap_nt_nonvar[i in G_nt_nonvar, t in T], 
        GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # Capacity constraints (variable generation)

    @constraint(UC, Cap_var[i in G_var, t in T], 
        GEN[i,t] <= gen_var[(gen_var.hour .== t) .& (gen_var.r_id .== i), :cf][1] *
                    gen_df[gen_df.r_id .== i, :existing_cap_mw][1])
    
    # battery constraints 
    @constraint(UC, cStateOfChargeEnd[b in B],
        SOC[b, T[end]] == 0.5 * storage_df[storage_df.r_id .== b, :existing_cap_mw][1]
    )
    @constraint(UC, cStateOfChargeStart[b in B],
        SOC[b, T[1]] == 0.5 * storage_df[storage_df.r_id .== b, :existing_cap_mw][1] + # start charge
                (CHARGE[b, T[1]] * storage_df[storage_df.r_id .== b, :charge_efficiency][1] -
                DISCHARGE[b, T[1]] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1])
    )
    @constraint(UC, cStateOfCharge[b in B, t in T[2:end]],
        SOC[b, t] == SOC[b, t-1] +
                    (CHARGE[b, t] * storage_df[storage_df.r_id .== b, :charge_efficiency][1] -
                    DISCHARGE[b, t] / storage_df[storage_df.r_id .== b, :discharge_efficiency][1])
    )

    # Unit commitment constraints
    @constraint(UC, Startup[i in G_thermal, t in T],
        COMMIT[i, t] >= sum(START[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:up_time][1]):t)))

    @constraint(UC, Shutdown[i in G_thermal, t in T],
        1-COMMIT[i, t] >= sum(SHUT[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:down_time][1]):t)))
    
    @constraint(UC, CommitmentStatus[i in G_thermal, t in T_red],
        COMMIT[i, t+1] - COMMIT[i, t] == START[i, t+1] - SHUT[i, t+1]
    )
    # New auxiliary variable GENAUX for generation above the minimum output level
    # for committed thermal units (only created for thermal generators)
    @constraint(UC, AuxGen[i in G_thermal, t in T],
        GENAUX[i,t] == GEN[i,t] - COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    
    # Ramp equations for thermal generators (constraining GENAUX)
    @constraint(UC, RampUp_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t+1] - GENAUX[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] )

    @constraint(UC, RampDn_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t] - GENAUX[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] )

    # Ramp equations for non-thermal generators (constraining total generation GEN)
    @constraint(UC, RampUp_nonthermal[i in G_nonthermal, t in T_red], 
        GEN[i,t+1] - GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] )

    @constraint(UC, RampDn[i in G, t in T_red], 
        GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] )
    
                                 
    # Reserve equations
    # (1) Reserves limited by committed capacity of generator
    @constraint(UC, ResUpCap[i in G_thermal, t in T],
        RESUP[i,t] <= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1]
                        - GEN[i,t])

    @constraint(UC, ResDnCap[i in G_thermal, t in T],
        RESDN[i,t] <= GEN[i,t] - 
                        COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])

    # (2) Reserves limited by ramp rates
    @constraint(UC, ResUpRamp[i in G_thermal, t in T],
        RESUP[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                          gen_df[gen_df.r_id .== i,:ramp_up_percentage][1])

    @constraint(UC, ResDnRamp[i in G_thermal, t in T],
        RESDN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                          gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])

    # Compute reserve requirements
    ResReqUp = Dict(T[i] => 300  
              + 0.05 * loads[loads.hour .== T[i],:demand][1]
                    for i in 1:length(T))
    ResReqDn = Dict(T[i] => 0
              + 0.05 * loads[loads.hour .== T[i],:demand][1]                    
                    for i in 1:length(T))
                        
    # (3) Overall reserve requirements
    @constraint(UC, ResUpRequirement[t in T],
        sum(RESUP[i,t] for i in G_thermal) >= 
                ResReqUp[t])

    @constraint(UC, ResDnRequirement[t in T],
        sum(RESDN[i,t] for i in G_thermal) >= 
                ResReqDn[t])
    
    # Solve statement (! indicates runs in place)
    optimize!(UC)

    # Generation solution
    gen = value_to_df_2dim(value.(GEN))

    # Commitment status solution
    commit = value_to_df_2dim(value.(COMMIT))

    # Calculate curtailment
    curtail = DataFrame(r_id = Int[], hour = Int[], curt = Float64[])
    for i in G_var
        for t in T
            cf = gen_var[(gen_var.hour .== t) .& (gen_var.r_id .== i), :cf][1]
            existing_cap = gen_df[gen_df.r_id .== i, :existing_cap_mw][1]
            generation = value(GEN[i,t])
            curt_val = cf * existing_cap - generation
            push!(curtail, (r_id = i, hour = t, curt = curt_val))
        end
    end

    # Calculate total CO2 emissions
    total_co2 = sum(gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    value(GEN[i,t]) for i in G_nonvar for t in T) +
                sum(gen_df[gen_df.r_id .== i,:start_fuel_mmbtu_per_mw][1] * 
                    gen_df[gen_df.r_id .== i,:co2_content_tons_per_mmbtu][1] * 
                    gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                    value(START[i,t]) for i in G_thermal for t in T)

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

    # Return the solution and objective as named tuple
    return (
        gen,
        commit,
        curtail,
        total_co2,
        charge_df,
        discharge_df,
        flows,
        reserves = value_to_df_2dim(value.(RESUP)),
        cost = objective_value(UC),
        status = termination_status(UC)
    )

end