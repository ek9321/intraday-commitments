using Plots, GraphRecipes, Graphs, Colors, Statistics, StatsBase

function get_zone_coordinates()
    return Dict(
        1 => (-122.4, 37.8),   # NorCal (San Francisco)
        2 => (-120.2, 34.1),   # SoCal (Los Angeles)
        3 => (-117.2, 32.7),   # SD_IID (San Diego)
        12 => (-104.9, 39.7),  # CO (Denver)
        8 => (-112.1, 33.4),   # AZ (Phoenix)
        9 => (-106.6, 35.1),   # NM (Albuquerque)
        5 => (-116.2, 43.6),   # ID (Boise)
        11 => (-111.9, 40.8),  # UT (Salt Lake City)
        7 => (-119.8, 39.5),   # NV (Reno)
        10 => (-106.3, 42.8),  # WY (Casper)
        4 => (-122.3, 47.6),   # PNW (Seattle)
        6 => (-108.5, 45.8)    # MT (Billings)
    )
end

"""
Create network graph with power flows
"""
function plot_power_network(
    flows_df,              # DataFrame with columns: line, hour, flow
    network,               # Network data with line definitions
    gen_df,                # Generator data with zone info
    zone_dict;             # Zone ID -> Zone name mapping
    hour=nothing,          # Specific hour to plot (or nothing for average)
    hours=nothing,         # All hours in dataset (for animation)
    title_str="Power Network Flows",
    node_size_by_gen=true, # Size nodes by generation capacity
    save_path="network_flows.png",
    arrow_scale=1.0
)
    
    coords = get_zone_coordinates()
    zones = sort(unique(gen_df.zone))
    n_zones = length(zones)
    
    # Filter flows for specific hour or average
    if hour !== nothing
        flows_hour = flows_df[flows_df.hour .== hour, :]
        if hours !== nothing
            hour = hour - minimum(hours)  # Normalize hour for title
        end
        title_str = "$title_str - Hour $hour"
    else
        # Average over all hours
        flows_hour = combine(groupby(flows_df, :line), :flow => mean => :flow)
        title_str = "$title_str - Average"
    end
    
    # Create adjacency matrix for graph
    adj_matrix = zeros(n_zones, n_zones)
    flow_matrix = zeros(n_zones, n_zones)
    
    for row in eachrow(flows_hour)
        line_idx = row.line
        flow_val = row.flow
        
        # Find which zones this line connects
        for i in zones
            for j in zones
                if i != j
                    # Check if line connects zones i and j
                    coef_i = network[line_idx, string(i)]
                    coef_j = network[line_idx, string(j)]
                    
                    if coef_i != 0 && coef_j != 0
                        adj_matrix[i, j] = 1
                        # Store flow (positive = i->j, negative = j->i)
                        if coef_i > 0  # Flow goes into zone i
                            flow_matrix[j, i] = abs(flow_val)
                        else
                            flow_matrix[i, j] = abs(flow_val)
                        end
                    end
                end
            end
        end
    end
    
    # Node sizes based on generation capacity
    if node_size_by_gen
        node_sizes = [
            sum(gen_df[gen_df.zone .== z, :existing_cap_mw]) / 1000  # Scale to GW
            for z in zones
        ]
        node_sizes = node_sizes ./ maximum(node_sizes) .* 20 .+ 10  # Scale for visualization
    else
        node_sizes = fill(10, n_zones)
    end
    
    # Node colors based on generation type dominance
    node_colors = []
    for z in zones
        gens_in_zone = gen_df[gen_df.zone .== z, :]
        ces_cap = sum(gens_in_zone[gens_in_zone.ces .== 1, :existing_cap_mw])
        total_cap = sum(gens_in_zone.existing_cap_mw)
        ces_fraction = ces_cap / total_cap
        
        # Color gradient: red (fossil) -> yellow (mixed) -> green (renewable)
        push!(node_colors, RGBA(1-ces_fraction, ces_fraction, 0, 0.8))
    end
    
    # Get node positions from coordinates
    node_x = [coords[z][1] for z in zones]
    node_y = [coords[z][2] for z in zones]
    
    # Create the plot
    p = plot(
        size=(1400, 1000),
        title=title_str,
        titlefontsize=16,
        legend=:topright,
        showaxis=false,
        grid=false,
        aspect_ratio=:equal
    )
    
    # Draw edges with arrows for each connection
    for row in eachrow(flows_hour)
        line_idx = row.line
        flow_val = abs(row.flow)
        
        if flow_val > 0
            # Find zones connected by this line
            for i in zones
                for j in zones
                    if i < j
                        coef_i = network[line_idx, string(i)]
                        coef_j = network[line_idx, string(j)]
                        
                        if coef_i != 0 && coef_j != 0
                            # Determine direction based on coefficients
                            if coef_i > 0
                                from_zone, to_zone = j, i
                            else
                                from_zone, to_zone = i, j
                            end
                            
                            # Edge width and color based on flow
                            max_flow = maximum(flows_hour.flow)
                            edge_width = (flow_val / max_flow) * 8 + 1.5
                            edge_alpha = min(1.0, (flow_val / max_flow) * 0.8 + 0.3)
                            
                            # Get coordinates
                            x1, y1 = node_x[findfirst(==(from_zone), zones)], node_y[findfirst(==(from_zone), zones)]
                            x2, y2 = node_x[findfirst(==(to_zone), zones)], node_y[findfirst(==(to_zone), zones)]
                            
                            # Calculate arrow direction and shorten line
                            dx = x2 - x1
                            dy = y2 - y1
                            dist = sqrt(dx^2 + dy^2)
                            node_radius = 0.5
                            
                            if dist > 2 * node_radius
                                scale = (dist - node_radius) / dist
                                x2_adj = x1 + dx * scale
                                y2_adj = y1 + dy * scale
                                
                                scale_start = node_radius / dist
                                x1_adj = x1 + dx * scale_start
                                y1_adj = y1 + dy * scale_start
                            else
                                x1_adj, y1_adj = x1, y1
                                x2_adj, y2_adj = x2, y2
                            end
                            
                            # Draw arrow with explicit arrow parameters
                            plot!(p, 
                                [x1_adj, x2_adj], 
                                [y1_adj, y2_adj],
                                arrow=(:closed, :head, 2.0, 2.0),  # Larger, more visible arrows
                                linewidth=edge_width,
                                color=RGBA(0.2, 0.4, 0.8, edge_alpha),
                                label=""
                            )
                            
                            # Add flow label at midpoint
                            mid_x = (x1 + x2) / 2
                            mid_y = (y1 + y2) / 2
                            annotate!(p, mid_x, mid_y, 
                                     text("$(round(flow_val, digits=0))", 8, :blue))
                        end 
                    end
                end
            end
        end
    end
    
    
    # Draw nodes
    scatter!(p,
        node_x,
        node_y,
        markersize=node_sizes,
        markercolor=node_colors,
        markerstrokewidth=2,
        markerstrokecolor=:black,
        label="",
        series_annotations=text.([zone_dict[z] for z in zones], 9, :black, :center)
    )
        
    # Add legend
    scatter!(p, [], [], 
            markersize=15, 
            markercolor=RGB(1, 0, 0),
            label="Fossil-heavy")
    scatter!(p, [], [],
            markersize=15,
            markercolor=RGB(0.5, 0.5, 0),
            label="Mixed")
    scatter!(p, [], [],
            markersize=15,
            markercolor=RGB(0, 1, 0),
            label="Renewable-heavy")

    if save_path !== nothing
        savefig(p, save_path)
    end

    return p
end

"""
Create animated network flows over time
"""
function animate_power_network(
    flows_df,
    network,
    gen_df,
    zone_dict;
    hours=nothing,  # Subset of hours to animate
    fps=2,
    save_path="network_animation.gif"
    )
    
    if hours === nothing
        hours = sort(unique(flows_df.hour))
    end

    # hours = hours .- minimum(hours)  # Normalize hours for animation
    anim = @animate for hour in hours
        plot_power_network(
            flows_df,
            network,
            gen_df,
            zone_dict;
            hour=hour,
            hours=hours,
            title_str="Power Network Flows",
            save_path=nothing
        )
    end
    
    gif(anim, save_path, fps=fps)
    println("Animation saved to $save_path")
end

"""
Create a heatmap showing flows on all lines over time
"""
function plot_flow_heatmap(
    flows_df,
    network,
    line_dict;  # Line ID -> Line name mapping
    save_path="flow_heatmap.png"
)
    
    hours = sort(unique(flows_df.hour))
    lines = sort(unique(flows_df.line))
    
    # Create matrix: lines x hours
    flow_matrix = zeros(length(lines), length(hours))
    
    for (i, line) in enumerate(lines)
        for (j, hour) in enumerate(hours)
            flow_val = - flows_df[(flows_df.line .== line) .& (flows_df.hour .== hour), :flow]
            if !isempty(flow_val)
                flow_matrix[i, j] = flow_val[1]
            end
        end
    end
    
    # Get line names
    line_names = [haskey(line_dict, l) ? line_dict[l] : "Line $l" for l in lines]
    
    p = heatmap(
        hours,
        1:length(lines),
        flow_matrix,
        xlabel="Hour",
        ylabel="Transmission Line",
        title="Power Flows Across Network",
        color=:RdBu,
        clims=(-maximum(abs.(flow_matrix)), maximum(abs.(flow_matrix))),
        yticks=(1:length(lines), line_names),
        size=(1200, 800),
        colorbar_title="Flow (MW)"
    )
    
    savefig(p, save_path)
    return p
end

"""
Compare flows before and after RUC
"""
function compare_network_flows(
    flows_uc,
    flows_ruc,
    network,
    gen_df,
    zone_dict;
    hour=nothing,
    save_path="flow_comparison.png"
)
    
    if hour !== nothing
        flows_uc_h = flows_uc[flows_uc.hour .== hour, :]
        flows_ruc_h = flows_ruc[flows_ruc.hour .== hour, :]
        title_suffix = " - Hour $hour"
    else
        flows_uc_h = combine(groupby(flows_uc, :line), :flow => mean => :flow)
        flows_ruc_h = combine(groupby(flows_ruc, :line), :flow => mean => :flow)
        title_suffix = " - Average"
    end
    
    p1 = plot_power_network(
        flows_uc,
        network,
        gen_df,
        zone_dict;
        hour=hour,
        title_str="Before RUC" * title_suffix,
        save_path=nothing,
        arrow_scale=0.2
    )
    
    p2 = plot_power_network(
        flows_ruc,
        network,
        gen_df,
        zone_dict;
        hour=hour,
        title_str="After RUC" * title_suffix,
        save_path=nothing,
        arrow_scale=0.2
    )
    
    # Calculate flow changes
    flows_change = innerjoin(
        flows_uc_h,
        flows_ruc_h,
        on=:line,
        makeunique=true
    )
    flows_change.flow_change = flows_change.flow_1 .- flows_change.flow
    
    # Plot just the changes
    coords = get_zone_coordinates()
    zones = sort(unique(gen_df.zone))
    node_x = [coords[z][1] for z in zones]
    node_y = [coords[z][2] for z in zones]
    
    p3 = plot(
        size=(1400, 1000),
        title="Flow Changes (RUC - UC)" * title_suffix,
        showaxis=false,
        grid=false,
        aspect_ratio=:equal
    )
    
    # Draw nodes
    scatter!(p3,
        node_x,
        node_y,
        markersize=15,
        markercolor=:lightgray,
        markerstrokewidth=2,
        markerstrokecolor=:black,
        label="",
        series_annotations=text.([zone_dict[z] for z in zones], 10, :black, :center)
    )
    
    # Draw flow changes
    for row in eachrow(flows_change)
        change = row.flow_change
        if abs(change) > 10  # Only show significant changes
            line_idx = row.line
            
            # Find zones connected by this line
            for i in zones
                for j in zones
                    if i < j
                        coef_i = network[line_idx, string(i)]
                        coef_j = network[line_idx, string(j)]
                        
                        if coef_i != 0 && coef_j != 0
                            # Determine direction
                            if change > 0
                                color = :green
                                arrow_start, arrow_end = i, j
                            else
                                color = :red
                                arrow_start, arrow_end = j, i
                            end
                            
                            plot!(p3,
                                [node_x[arrow_start], node_x[arrow_end]],
                                [node_y[arrow_start], node_y[arrow_end]],
                                arrow=true,
                                linewidth=min(abs(change)/50, 10),
                                color=color,
                                alpha=0.6,
                                label=""
                            )
                            
                            mid_x = (node_x[i] + node_x[j]) / 2
                            mid_y = (node_y[i] + node_y[j]) / 2
                            annotate!(p3, mid_x, mid_y,
                                     text("$(round(change, digits=0))", 8, color))
                        end
                    end
                end
            end
        end
    end
    
    # Combine all three plots
    plot(p1, p2, p3, layout=(3,1), size=(1400, 2400))
    savefig(save_path)
end

"""
Create congestion analysis plot
"""
function plot_congestion_analysis(
    flows_df,
    network;
    save_path="congestion_analysis.png"
)
    
    lines = sort(unique(flows_df.line))
    hours = sort(unique(flows_df.hour))
    
    # Calculate congestion metrics
    congestion_data = DataFrame(
        line=Int[],
        max_capacity=Float64[],
        avg_flow=Float64[],
        peak_flow=Float64[],
        utilization=Float64[],
        hours_congested=Int[]
    )
    
    for line in lines
        line_flows = flows_df[flows_df.line .== line, :flow]
        max_cap = network.line_max_flow_mw[line]
        
        avg_flow = mean(abs.(line_flows))
        peak_flow = maximum(abs.(line_flows))
        utilization = peak_flow / max_cap
        hours_congested = sum(abs.(line_flows) .> 0.9 * max_cap)
        
        push!(congestion_data, (
            line=line,
            max_capacity=max_cap,
            avg_flow=avg_flow,
            peak_flow=peak_flow,
            utilization=utilization,
            hours_congested=hours_congested
        ))
    end
    
    # Sort by hours congested
    sort!(congestion_data, :hours_congested, rev=true)
    
    # Plot top 10 most frequently congested lines
    top_n = min(10, nrow(congestion_data))
    
    p = bar(
        1:top_n,
        congestion_data[1:top_n, :hours_congested],
        xlabel="Transmission Line",
        ylabel="Hours Above 90% Capacity",
        title="Most Frequently Congested Transmission Lines (>90% Capacity)",
        legend=false,
        color=ifelse.(congestion_data[1:top_n, :hours_congested] .> 0, :red, :green),
        xticks=(1:top_n, [haskey(line_dict, congestion_data[i, :line]) ? line_dict[congestion_data[i, :line]] : "Line $(congestion_data[i, :line])" for i in 1:top_n]),
        xrotation=45,
        size=(1200, 600)
    )
    
    savefig(p, save_path)
    return congestion_data
end