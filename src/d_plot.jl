# Load necessary libraries
using Plots
using Plots.Measures
using DelimitedFiles
using Printf

# Select the plotting backend
gr()

# --- Function to Plot the Data ---
function plot_all_data_subplots(filepath::String)
    
    # Load data from the file, ignoring comment lines
    data = readdlm(filepath, comments=true)
    
    if isempty(data)
        println("File is empty: $filepath. Cannot plot.")
        return
    end

    # Extract all relevant columns
    d       = data[:, 1]
    m       = data[:, 2]
    q       = data[:, 3]
    Q       = data[:, 5]
    q_bar   = data[:, 6]
    λ_asymm = data[:, 7]
    λ_symm = data[:, 8]
    λ_unc = data[:, 9]
    
    # --- PLOT 1 CREATION: Order Parameters (left subplot) ---
    plot1 = plot(d, m, 
        label="m", 
        legend=:topright, 
        lw=2.5, 
        color=:blue,
        ylabel="Order Parameters")
        
    plot!(plot1, d, q, label="q", lw=2, color=:purple)
    plot!(plot1, d, q_bar, label="q̄", lw=2, color=:darkgreen, linestyle=:dash)
    plot!(plot1, d, Q, label="Q", lw=2, color=:darkorange, linestyle=:dot)
    xlabel!(plot1, "Dilution d") # Add x-label only to the plots on the bottom row if needed

    # --- PLOT 2 CREATION: Eigenvalues (right subplot) ---
    plot2 = plot(d, λ_asymm, 
        label="λ_asymm", 
        legend=:topright, 
        lw=2, 
        color=:red,
        ylabel="Eigenvalue λ")
        
    plot!(plot2, d, λ_symm, label="λ_symm", lw=2, color=:brown)
    plot!(plot2, d, λ_unc, label="λ_unc", lw=2, color=:cyan)
    
    # Add the stability line, but remove it from the legend with label=""
    hline!(plot2, [0], linestyle=:dash, color=:black, label="")
    xlabel!(plot2, "Dilution d")

    
    # --- Combine the two subplots into a single figure ---
    title_str = basename(filepath)
    main_title = "Analysis vs Dilution - $(splitext(title_str)[1])"
    
    # layout=(1, 2) arranges the plots in 1 row and 2 columns
    final_plot = plot(plot1, plot2, layout=(1, 2), size=(1200, 500), plot_title=main_title, bottom_margin=5mm, left_margin=5mm)
    
    # Display and save the plot
    display(final_plot)

    output_folder = "plots_scan_d_broad"
    !isdir(output_folder) && mkdir(output_folder)
    output_filename = joinpath(output_folder, "plot_subplots_side_$(splitext(title_str)[1]).png")
    savefig(final_plot, output_filename)
    println("Plot saved to: $output_filename")
end

# --- Execution ---
# This part remains the same. It finds and calls the plotting function.
s = 3
T_fixed = 0.001
alpha_range = [0.0]#0.0001:0.0001:0.02 #0.001:0.001:0.02

for alpha in alpha_range

    T_formatted = @sprintf "%.3f" T_fixed
    alpha_formatted = @sprintf "%.4f" alpha
    filename = "scan_d_broad/scan_d_s=$(s)_T=$(T_formatted)_alpha=$(alpha_formatted).txt"

    if isfile(filename)
        plot_all_data_subplots(filename)
    else
        println("File not found: $filename. Run the simulation script first.")
    end

end
