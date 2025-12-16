# Load necessary libraries
using Plots
using DelimitedFiles
using Printf

# Select the plotting backend
gr()
# pythonplot()

# --- PLOTTING FUNCTION DEFINITION ---
# This function creates a single plot and saves it to a specific subfolder.

function create_and_save_plot(s::Int, d::Float64, T_plot::Float64, output_basedir::String)
    
    println("Processing: s=$s, d=$d, T=$T_plot")

    # --- Data Loading and Preparation ---
    datapath = "result" # Path to the results folder
    # Format 'd' with zero-padding for correct alphabetical sorting
    d_formatted = @sprintf "%.3f" d
    filename = joinpath(datapath, "mix_s=$(s)_d=$(d_formatted).txt")

    if !isfile(filename)
        println("--> File not found: $filename. Skipping.")
        return # Exit the function if the file doesn't exist
    end

    data = readdlm(filename, comments=true)

    # Check if the file is empty
    if isempty(data)
        println("--> File is empty: $filename. Skipping.")
        return
    end

    # Extract columns based on the header in your data file
    α_full = data[:, 1]
    T_full = data[:, 2]
    m_full = data[:, 3]
    λ_asymm_full = data[:, 8]
    λ_symm_full = data[:, 9]
    λ_unc_full = data[:, 10]

    # Filter data for the chosen temperature.
    # Due to floating-point inaccuracies, it's safer to use isapprox (≈)
    indices = findall(t -> isapprox(t, T_plot, atol=1e-4), T_full)

    if isempty(indices)
        println("--> No data found for T=$T_plot in file $filename. Skipping.")
        return
    end

    α = α_full[indices]
    m = m_full[indices]
    λ_asymm = λ_asymm_full[indices]
    λ_symm = λ_symm_full[indices]
    λ_unc = λ_unc_full[indices]

    # --- Combined Plot Creation ---
    plot_left = plot(α, m, label="m", xlabel="Load α", ylabel="Magnetization m", color=:blue, legend=:bottomleft, lw=3)
    
    plot_right = twinx(plot_left)
    
    plot!(plot_right, α, λ_asymm, label="λ_asymm", ylabel="Eigenvalue λ", color=:red, lw=2)
    plot!(plot_right, α, λ_symm, label="λ_symm", color=:green, lw=2)
    plot!(plot_right, α, λ_unc, label="λ_unc", color=:orange, lw=2)
    hline!(plot_right, [0], linestyle=:dash, color=:black, label="")
    title!("s=$s, d=$d, T=$T_plot")
    
    # --- Plot Saving ---
    # Create the subfolder for the current temperature if it doesn't exist
    temp_folder = joinpath(output_basedir, "T_$(T_plot)")
    !isdir(temp_folder) && mkdir(temp_folder)
    
    output_filename = joinpath(temp_folder, "plot_s=$(s)_d=$(d_formatted).png")
    savefig(plot_left, output_filename)
    println("--> Plot saved to: $output_filename")
end


# --- MAIN EXECUTION SCRIPT ---

# --- Experiment Parameters ---
s_val = 3
output_directory = "plots_by_temperature" # Main output folder for all plots

# List of temperatures to create folders and plots for
T_values = [0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001]

# List of 'd' values to analyze
d_values = 0.0:0.025:1.0

# Create the main output directory if it doesn't exist
!isdir(output_directory) && mkdir(output_directory)


# --- Execution Loop ---
# Double loop: one over temperature, one over dilution
for T_plot in T_values
    println("\n--- Generating all plots for T = $T_plot ---")
    for d_val in d_values
        # Call the plotting function with the current scalar values
        create_and_save_plot(s_val, d_val, T_plot, output_directory)
    end
end

println("\nAll plots generated!")

