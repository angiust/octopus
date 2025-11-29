# Carica tutti i tuoi moduli. Assicurati che i path siano corretti.
# include("hopfield.jl")
# include("mixtures.jl")
include("spinodal.jl")

# --- Definizione dei Parametri per la Scansione ---

# 1. Parametri della miscela
s = 3 # Miscela di 3 pattern

# 2. Range per la diluizione 'd'
#    Va da 0.0 a 1.0 con passo 0.1
d_range = 0.0:0.025:1.0

# 3. Range per la temperatura 'T'
#    Va da 0.05 a 0.5 con passo 0.05 (ho usato 0.05 invece di 0.005 per un test più veloce)
#    Se vuoi un passo di 0.005, usa: T_range = 0.05:0.005:0.5
T_range = 0.01 #0.005:0.005:0.2
beta_range = 1 ./ T_range

# 4. Range per il carico 'α'
#    (Usiamo gli stessi valori del tuo test precedente)
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
