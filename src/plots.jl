function plot_generation_by_zone(solution, T_period, gen_df, zone_dict, title = "Generation by Technology Across Zones")
    # Plot generation by technology for each zone in a 4x3 grid
    # Join generator zone info
    sol_gen_zone = innerjoin(solution.gen, gen_df[:, [:r_id, :technology, :zone]], on = :r_id)

    # Group by technology, zone, and hour, then sum generation
    sol_gen_zone = combine(groupby(sol_gen_zone, [:technology, :zone, :hour]), :gen => sum => :gen_sum)

    # Add zone names
    sol_gen_zone = innerjoin(sol_gen_zone, DataFrame(zone = collect(keys(zone_dict)), zone_name = collect(values(zone_dict))), on = :zone)

    # Rescale hours
    sol_gen_zone.hour = sol_gen_zone.hour .- T_period[1]

    # Create the grid plot with consistent colors
    sol_gen_zone |>
    @vlplot(
        columns = 3,
        facet = {field = "zone_name", type = "nominal"},
        spec = {
            mark = :area,
            width = 200,
            height = 150,
            encoding = {
                x = {field = :hour, type = "quantitative"},
                y = {field = :gen_sum, type = "quantitative", stack = :zero, title = "Generation (MW)"},
                color = {field = "technology", type = "nominal", scale = {scheme = "category10"}}
            }
        },
        title = title,
        resolve = {scale = {color = "independent"}}
    ) |> display

    return sol_gen_zone
end

function plot_ruc_generation_by_zone(solution, ruc_solution, initial_time, recommit_time, gen_df, zone_dict, title)
    # Plot generation by technology for each zone in a 4x3 grid for RUC solution
    # Combined with original solution for first 12 hours

    # Combine original and RUC solutions by zone
    combined_gen_zone = DataFrame()

    # Original solution up to recommit_time (exclusive)
    original_cutoff_zone = solution.gen[solution.gen.hour .< (initial_time + recommit_time), :]
    original_cutoff_zone = innerjoin(original_cutoff_zone, gen_df[:, [:r_id, :technology, :zone]], on = :r_id)
    original_cutoff_zone = combine(groupby(original_cutoff_zone, [:technology, :zone, :hour]), :gen => sum => :gen_sum)

    # RUC solution from recommit_time onwards (inclusive)
    ruc_gen_zone = innerjoin(ruc_solution.gen, gen_df[:, [:r_id, :technology, :zone]], on = :r_id)
    ruc_gen_zone = combine(groupby(ruc_gen_zone, [:technology, :zone, :hour]), :gen => sum => :gen_sum)

    # Combine both parts
    append!(combined_gen_zone, original_cutoff_zone)
    append!(combined_gen_zone, ruc_gen_zone)

    # Add zone names
    combined_gen_zone = innerjoin(combined_gen_zone, DataFrame(zone = collect(keys(zone_dict)), zone_name = collect(values(zone_dict))), on = :zone)

    # Rescale hours
    combined_gen_zone.hour = combined_gen_zone.hour .- initial_time

    # Create the grid plot with consistent colors
    combined_gen_zone |>
    @vlplot(
        columns = 3,
        facet = {field = "zone_name", type = "nominal"},
        spec = {
            mark = :area,
            width = 200,
            height = 150,
            encoding = {
                x = {field = :hour, type = "quantitative"},
                y = {field = :gen_sum, type = "quantitative", stack = :zero, title = "Generation (MW)"},
                color = {field = "technology", type = "nominal", scale = {scheme = "category10"}}
            }
        },
        title = "Generation by Technology with RUC Interruption at Hour $(recommit_time)",
        resolve = {scale = {color = "independent"}}
    ) |> display

    return combined_gen_zone
end

function plot_total_generation(solution, T_period, gen_df, zone_dict, title)
    sol_gen = innerjoin(solution.gen, 
                    gen_df[!, [:r_id, :technology]], 
                    on = :r_id)
    sol_gen = combine(groupby(sol_gen, [:technology, :hour]), 
                :gen => sum)

    # Curtailment
    curtail = combine(groupby(solution.curtail, [:hour]),
                :curt => sum)
    curtail[!, :technology] .= "_curtailment"
    rename!(curtail, :curt_sum => :gen_sum)
    append!(sol_gen, curtail[:,[:technology, :hour, :gen_sum]])

    # Add storage discharge
    for (storagetype, storagelabel, df) in [
        ("discharge", "_storage_discharge", solution.discharge_df),
        ("charge", "_storage_charge", solution.charge_df)
        ]
        storage_plot = combine(groupby(df, :hour), :value => sum)
        storage_plot[!, :technology] .= storagelabel
        rename!(storage_plot, :value_sum => :gen_sum)
        append!(sol_gen, storage_plot[:, [:technology, :hour, :gen_sum]])
    end

    # Rescale hours
    sol_gen.hour = sol_gen.hour .- T_period[1]
    # Add plot title
    title = "Generation by Technology (No RUC)"
    sol_gen |>
    @vlplot(:area, 
        x=:hour, y={:gen_sum, stack=:zero}, 
        color={"technology:n", scale={scheme="category10"}}
    ) |> display

    return sol_gen
end

function plot_total_ruc_generation(solution, ruc_solution, initial_time, recommit_time, gen_df, zone_dict, title)
    # Combine original and RUC solutions with interruption at recommit_time
    combined_gen = DataFrame()

    # Original solution up to recommit_time (exclusive)
    original_cutoff = solution.gen[solution.gen.hour .< (initial_time + recommit_time), :]
    original_cutoff = innerjoin(original_cutoff, gen_df[!, [:r_id, :technology]], on = :r_id)
    original_cutoff = combine(groupby(original_cutoff, [:technology, :hour]), :gen => sum)

    # RUC solution from recommit_time onwards (inclusive)
    ruc_gen = innerjoin(ruc_solution.gen, gen_df[!, [:r_id, :technology]], on = :r_id)
    ruc_gen = combine(groupby(ruc_gen, [:technology, :hour]), :gen => sum)

    # Combine both parts
    append!(combined_gen, original_cutoff)
    append!(combined_gen, ruc_gen)

    # Add curtailment from both solutions
    original_curtail = solution.curtail[solution.curtail.hour .< (initial_time + recommit_time), :]
    original_curtail = combine(groupby(original_curtail, [:hour]), :curt => sum)
    original_curtail[!, :technology] .= "_curtailment"
    rename!(original_curtail, :curt_sum => :gen_sum)

    ruc_curtail = combine(groupby(ruc_solution.curtail, [:hour]), :curt => sum)
    ruc_curtail[!, :technology] .= "_curtailment"
    rename!(ruc_curtail, :curt_sum => :gen_sum)

    append!(combined_gen, original_curtail[:,[:technology, :hour, :gen_sum]])
    append!(combined_gen, ruc_curtail[:,[:technology, :hour, :gen_sum]])

    # Add storage for both solutions
    for (sol, time_filter) in [(solution, x -> x < (initial_time + recommit_time)), 
                            (ruc_solution, x -> true)]
        for (storagetype, storagelabel, df) in [
            ("discharge", "_storage_discharge", sol.discharge_df),
            ("charge", "_storage_charge", sol.charge_df)
            ]
            storage_data = df[time_filter.(df.hour), :]
            storage_plot = combine(groupby(storage_data, :hour), :value => sum)
            storage_plot[!, :technology] .= storagelabel
            rename!(storage_plot, :value_sum => :gen_sum)
            append!(combined_gen, storage_plot[:, [:technology, :hour, :gen_sum]])
        end
    end

    # Rescale hours
    combined_gen.hour = combined_gen.hour .- initial_time

    combined_gen |>
    @vlplot(:area, 
        x=:hour, y={:gen_sum, stack=:zero}, 
        color={"technology:n", scale={scheme="category10"}},
        title="Generation with RUC Interruption at Hour $(recommit_time)"
        ) |> display

    return combined_gen
end

function plot_actual_generation(solution, T_period, gen_df, zone_dict, title)
    sol_gen = innerjoin(solution.gen, 
                    gen_df[!, [:r_id, :technology]], 
                    on = :r_id)
    sol_gen = combine(groupby(sol_gen, [:technology, :hour]), 
                :gen => sum)

    # Curtailment
    curtail = combine(groupby(solution.curtail, [:hour]),
                :curt => sum)
    curtail[!, :technology] .= "_curtailment"
    rename!(curtail, :curt_sum => :gen_sum)
    append!(sol_gen, curtail[:,[:technology, :hour, :gen_sum]])

    # Load shedding
    loadshed = combine(groupby(solution.loadshed, [:hour]),
                :gen => sum)  # assuming column is :gen in loadshed df
    loadshed[!, :technology] .= "_load_shedding"
    rename!(loadshed, :gen_sum => :gen_sum)
    append!(sol_gen, loadshed[:,[:technology, :hour, :gen_sum]])

    # Add storage discharge
    for (storagetype, storagelabel, df) in [
        ("discharge", "_storage_discharge", solution.discharge_df),
        ("charge", "_storage_charge", solution.charge_df)
        ]
        storage_plot = combine(groupby(df, :hour), :value => sum)
        storage_plot[!, :technology] .= storagelabel
        rename!(storage_plot, :value_sum => :gen_sum)
        append!(sol_gen, storage_plot[:, [:technology, :hour, :gen_sum]])
    end

    # Rescale hours
    sol_gen.hour = sol_gen.hour .- T_period[1]
    # Add plot title
    title = "Generation by Technology (No RUC)"
    sol_gen |>
    @vlplot(:area, 
        x=:hour, y={:gen_sum, stack=:zero}, 
        color={"technology:n", scale={scheme="category10"}}
    ) |> display

    return sol_gen
end

function plot_actual_generation_by_zone(solution, T_period, gen_df, zone_dict, title = "Generation by Technology Across Zones (Actual Conditions)")
    # Join generator zone info
    sol_gen_zone = innerjoin(solution.gen, gen_df[:, [:r_id, :technology, :zone]], on = :r_id)
    
    # Group by technology, zone, and hour, then sum generation
    sol_gen_zone = combine(groupby(sol_gen_zone, [:technology, :zone, :hour]), :gen => sum => :gen_sum)
    
    # Add curtailment by zone
    curtail_zone = innerjoin(solution.curtail, gen_df[:, [:r_id, :zone]], on = :r_id)
    curtail_zone = combine(groupby(curtail_zone, [:zone, :hour]), :curt => sum => :gen_sum)
    curtail_zone[!, :technology] .= "_curtailment"
    append!(sol_gen_zone, curtail_zone)
    
    # Add load shedding by zone
    loadshed_zone = copy(solution.loadshed)
    # loadshed has :r_id column which represents zone, rename to :zone
    rename!(loadshed_zone, :r_id => :zone)
    loadshed_zone = combine(groupby(loadshed_zone, [:zone, :hour]), :gen => sum => :gen_sum)
    loadshed_zone[!, :technology] .= "_load_shedding"
    append!(sol_gen_zone, loadshed_zone[:, [:technology, :zone, :hour, :gen_sum]])
    
    # Add zone names
    sol_gen_zone = innerjoin(sol_gen_zone, DataFrame(zone = collect(keys(zone_dict)), zone_name = collect(values(zone_dict))), on = :zone)
    
    # Rescale hours
    sol_gen_zone.hour = sol_gen_zone.hour .- T_period[1]
    
    # Create the grid plot
    sol_gen_zone |>
    @vlplot(
        columns = 3,
        facet = {field = "zone_name", type = "nominal"},
        spec = {
            mark = :area,
            width = 200,
            height = 150,
            encoding = {
                x = {field = :hour, type = "quantitative"},
                y = {field = :gen_sum, type = "quantitative", stack = :zero, title = "Generation (MW)"},
                color = {field = "technology", type = "nominal", scale = {scheme = "category10"}}
            }
        },
        title = title,
        resolve = {scale = {color = "independent"}}
    ) |> display
    
    return sol_gen_zone
end