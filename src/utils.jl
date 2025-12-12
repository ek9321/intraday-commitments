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