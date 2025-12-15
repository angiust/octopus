# include("hopfield.jl")
# include("mixtures.jl")
include("spinodal.jl")

# --- Definizione dei Parametri per la Scansione ---

# 1. Parametri della miscela
s = 3 # Miscela di 3 pattern

# 2. Range per la diluizione 'd'
d_range = 0.0:0.025:1.0

# 3. Range per la temperatura 'T'
T_range = [0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001] #0.005:0.005:0.2
beta_range = 1 ./ T_range

# 4. Range per il carico 'α'
alpha_range = 0.001:0.001:0.2


# --- Funzione Principale che Esegue Tutti gli Esperimenti ---

function run_all_scans()
    println("===== STARTING FULL PARAMETER SCAN =====")
    
    # Loop esterno sulla diluizione 'd'
    for d_val in d_range
        println("\n\n<<<<< RUNNING FOR DILUTION d = $d_val >>>>>\n")
        
        # Per ogni valore di d, chiama la tua funzione mix_span.
        # mix_span si occuperà dei loop interni su α e β.
        mix_span(
            collect(alpha_range), # Vettore di alpha
            collect(beta_range),  # Vettore di beta
            s,                    # Numero di pattern nella miscela
            d_val,                # Valore corrente di diluizione
            #resfile_prefix = "./results/data" # Prefisso per i file di output
        )
        
        println("\n<<<<< FINISHED FOR d = $d_val >>>>>")
    end
    
    println("\n===== ALL SCANS COMPLETED =====")
end

# --- Esecuzione ---

# Chiama la funzione principale per lanciare tutti i calcoli
run_all_scans()
