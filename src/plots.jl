# Carica le librerie necessarie
using Plots
using DelimitedFiles
using Glob # Utile per trovare i file

# Seleziona il backend per i plot
gr()
# pythonplot()

# --- Definizione della Funzione di Plotting ---
# Ho incapsulato il tuo codice in una funzione che prende
# i parametri s, d, e T come input.

function create_and_save_plot(s::Int, d::Float64, T_plot::Float64)
    
    println("Processing plot for s=$s, d=$d, T=$T_plot...")

    # --- Caricamento e Preparazione dei Dati ---
    datapath = "result" # Path alla cartella dei risultati
    filename = joinpath(datapath, "mix_s=$(s)_d=$(d).txt")

    if !isfile(filename)
        println("--> File not found: $filename. Skipping.")
        return # Esce dalla funzione se il file non esiste
    end

    data = readdlm(filename, comments=true)

    # Controlla se il file è vuoto
    if isempty(data)
        println("--> File is empty: $filename. Skipping.")
        return
    end

    # Estrai le colonne
    α_full = data[:, 1]
    T_full = data[:, 2]
    m_full = data[:, 3]
    λ_asymm_full = data[:, 8]
    λ_symm_full = data[:, 9]
    λ_unc_full = data[:, 10]

    # Filtra i dati per la temperatura scelta.
    # A causa di imprecisioni floating-point, è meglio usare isapprox (≈)
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

    # --- Creazione del Plot Combinato (il tuo codice è perfetto qui) ---
    plot_left = plot(α, m,
        label="m (Magnetization)",
        xlabel="Load α",
        ylabel="Magnetization m",
        color=:blue,
        legend=:bottomleft,
        lw=3)

    plot_right = twinx(plot_left)

    plot!(plot_right, α, λ_asymm, label="λ_asymm", ylabel="Eigenvalue λ", color=:red, lw=2)
    plot!(plot_right, α, λ_symm, label="λ_symm", color=:green, lw=2)
    plot!(plot_right, α, λ_unc, label="λ_unc", color=:orange, lw=2)
    hline!(plot_right, [0], linestyle=:dash, color=:black, label="")
    title!("m and Eigenvalues vs α (T=$T_plot, s=$s, d=$d)")
    
    # --- Salvataggio del Plot ---
    output_folder = "plots"
    !isdir(output_folder) && mkdir(output_folder) # Crea la cartella se non esiste
    
    output_filename = joinpath(output_folder, "plot_s=$(s)_d=$(d)_T=$(T_plot).png")
    savefig(plot_left, output_filename)
    println("--> Plot saved to: $output_filename")

end


# --- Esecuzione del Loop ---

# Parametri fissi per la serie di plot
s_val = 3
T_val = 0.01 # Plottiamo sempre la stessa temperatura T=0.2

# Lista di tutti i valori di 'd' per cui hai generato dati
d_values = collect(0.0:0.025:1.0) #[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]

# Loop che chiama la funzione di plotting per ogni valore di d
for d_val in d_values
    create_and_save_plot(s_val, d_val, T_val)
end

println("\nAll plots generated!")