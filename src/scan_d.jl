include("hopfield.jl")
include("mixtures.jl")
include("spinodal.jl")

using Printf


function run_scan_on_dilution()
    
    # parameters
    s = 3
    T_fixed = 0.001
    alpha_range = [0.0]#0.0001:0.0001:0.02
    d_range = 0.0:0.001:0.999

    for alpha_val in alpha_range

    	# folder and output filename
	output_dir = "scan_d_broad"
	!isdir(output_dir) && mkdir(output_dir)
	T_formatted = @sprintf "%.3f" T_fixed
	alpha_formatted = @sprintf "%.4f" alpha_val
	resfile = joinpath(output_dir, "scan_d_s=$(s)_T=$(T_formatted)_alpha=$(alpha_formatted).txt")

	println("===== STARTING SCAN ON DILUTION 'd' =====")
	println("Fixed parameters: s=$s, T=$T_fixed, α=$alpha_val")

	# header of output file
	open(resfile, "w") do f
	println(f,"# s=$s, T=$T_fixed, α=$alpha_val | 1:d 2:m 3:q 4:r 5:Q 6:q_bar 7:λ_asymm 8:λ_symm 9:λ_unc 10:converged(bool)")
	end

	# initialization
	pars = P.Params(ψ = 0.9, verb=0, ϵ=1e-7)
	op = P.OrderParams(1.0, 1, 0.5)

	# loop on dilution
	for d_val in d_range

	    # set the external parameters
	    ep = P.ExtParams(α=alpha_val, β=1/T_fixed, d=d_val)

	    # compute distribution of ζ for (s,d)
	    zp = calc_zeta_prob(s, d_val)

	    # convergence
	    ok, Δ = converge!(op, ep, pars, zp)

	    # compute quantities
	    Q = get_Q(op, ep, zp)
	    q_bar = get_q_bar(op, ep, zp)
	    λ_asymm = get_λ3(op, ep, zp, Q, q_bar)
	    λ_symm = get_λ1(op, ep, zp, Q, q_bar)
	    λ_unc = get_λ2(op, ep)

	    # save data
	    open(resfile, "a") do rf
	        println(rf,"$(ep.d) $(op.m) $(op.q) $(op.qh) $Q $q_bar $λ_asymm $λ_symm $λ_unc $ok")
	    end
	end

    end
    println("\n===== SCAN COMPLETED. Results saved =====")
end


run_scan_on_dilution()
