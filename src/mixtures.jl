using Parameters
using ExtractMacro

include("../common/common.jl")
include("./hopfield.jl")

#### struct per fare il valor medio sulla distribuzione di ζₛ ####

@with_kw mutable struct ZetaProb
    ζ = Int[]
    prob = Float64[]
    s::Int = 0
    # d::Float64 = 0.0 should i add it?
end

function calc_zeta_prob(s::Int, d::Float64)
    zp = ZetaProb()
    zp.s = s
    # zp.d = d it depends if i add it to the struct or not
    zp.ζ = collect(-s:s)

    zp.prob = map(x -> begin
        p = 0.0
        for k in 0:s # k = number of non-zero variables
            if (x + k) % 2 == 0 && abs(x) <= k
                prob_k_active = binomial(s, k) * (1-d)^k * d^(s-k)
                num_plus = (x + k) ÷ 2 # number of +1 in the k active variables
                prob_sum_is_x_given_k = binomial(k, num_plus) / 2.0^k
                p += prob_k_active * prob_sum_is_x_given_k
            end
        end
        p
    end, zp.ζ)

    return zp
end

#### free energy ####
function Gi_s(q, qh, m, α, β, s)
    α / 2 + 0.5 * s * m^2 + 0.5 * α * β * qh * (1 - q)
end

Gi_s(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = Gi_s(op.q, op.qh, op.m, ep.α, ep.β, zp.s)



#### ENERGETIC TERM ####
function Ge_s(qh, m, α, β, d, s, ζ, prob)
    # s = ζ[end] now i am passing s as a parameter of the function
    P.∫D(z -> begin
        result = 0.0
        for i in eachindex(ζ) #### valore medio sulla distribuzione di ζ
            term = √(α * qh) * z + m * ζ[i] 
            result += -log(2 * cosh(β * term)) * prob[i]
        end
        result / β # now i divided by beta and i don't need to multiply by 2 because the sum is over all the possible values of ζ
    end)
end

Ge_s(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = Ge_s(op.qh, op.m, ep.α, ep.β, ep.d, zp.s, zp.ζ, zp.prob)


function free_energy(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb)
    Gi_s(op, ep, zp) + Gs(op, ep) + Ge_s(op, ep, zp)
end


#### SADDLE POINT EQUATIONS p-MIXTURES RETRIEVAL ####

function update_sm(qh, m, α, β, s, ζ, prob) 
    # @assert length(prob) == length(ζ)
    # s = ζ[end]        
    P.∫D(z -> begin
        result = 0.0
        for i in eachindex(ζ) #### valore medio sulla distribuzione di ζ
            term = √(α * qh) * z + ζ[i] * m
            result += prob[i] * ζ[i] * tanh(β * term) /s # should i divide by s here or outside the integral?
        end
        result # now i don't need to multiply by 2 because the sum is over all the possible values of ζ
    end)
end

update_sm(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = update_sm(op.qh, op.m, ep.α, ep.β, zp.s, zp.ζ, zp.prob)

function update_sq(qh, m, α, β, ζ, prob) 
    # @assert length(prob) == length(ζ)
    # s = ζ[end]        
    P.∫D(z -> begin
        result = 0.0
        for i in eachindex(ζ) #### valore medio sulla distribuzione di ζ
            term = √(α * qh) * z + ζ[i] * m
            result += prob[i] * (tanh(β * term)^2)
        end
        result # now i don't need to multiply by 2 because the sum is over all the possible values of ζ
    end)
end
update_sq(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = update_sq(op.qh, op.m, ep.α, ep.β, zp.ζ, zp.prob)


#################  SADDLE POINT  ##################
# Right-hand-side
fsq(op, ep, zp) = update_sq(op, ep, zp)			   	   # q = fq (der: qh)
fsm(op, ep, zp) = update_sm(op, ep, zp)			   	   # m = fm (der: m)




function converge!(op::P.OrderParams, ep::P.ExtParams, pars::P.Params, zp::ZetaProb)
    @extract pars: maxiters verb ϵ ψ
    Δ = Inf
    ok = false   	
    

    for it = 1:maxiters
        Δ = 0.0
        ok = true
        verb > 1 && println("########## it=$it ##########")

        ########################################################################

        @update  op.qh   P.fqh       Δ ψ verb  op ep  #update qh
        @update  op.q    fsq       Δ ψ verb  op ep zp  #update q
        @update  op.m    fsm       Δ ψ verb  op ep zp  #update m
        
        
        ########################################################################

        verb > 1 && println(" Δ=$Δ\n")
        #(println(op); println(ep))
        #verb > 2 && it%5==0 && (println(ep);println(all_therm_func(op, ep));println(op))

        @assert isfinite(Δ)
        ok &= Δ < ϵ
        ok && break         # if ok==true, exit
    end
    ok, Δ
end


function calc_alpha_s(s::Int, d::Float64)

    zp = calc_zeta_prob(s, d)

    #### CALCOLO ARRAY β ####
    T_max = 0.50 #1. /sqrt(s) s = [3, 5, 7, 9]-mixtures, T = [0.45, 0.348, 0.29, 0.26] ---> usare questi valori
    T_min = 0.36
    δ_T = (T_max-T_min)/20 #divido per il numero di punti che voglio
    T = T_min:δ_T:T_max
    β = map(x -> 1. /x, T)

    #### CALCOLO ARRAY α ####
    α_max = 0.016 #dove finisce la linea spinoidale 3-mixtures. 0.138/s
    δ_α = α_max/1e3 #per le misture δ_α deve essere almeno ∝ 1e-6
    α = 1e-10:δ_α:α_max

    ep = P.ExtParams()
    ep.α = first(α)
    ep.β = first(β)
    ep.d = d

    op = P.OrderParams()
    op.m = 0.5#binomial(zp.s-1,convert(Int, (zp.s-1)/2))/(2^zp.s) viene riinizializzato dentro calc_alpha_s!

    pars = P.Params()

    resfile = "../results/data/mixtures_s=$(s)_d=$(d).txt"
    return calc_alpha_s!(op, ep, pars, zp, resfile, β, α)

end

function calc_alpha_s!(op::P.OrderParams, ep::P.ExtParams, pars::P.Params, zp::ZetaProb, resfile, β, α)

    !isfile(resfile) && open(resfile, "w") do f
    println(f,"# s = $(zp.s) d = $(ep.d) | 1:α_S 2:T")
    end

    m_old = 0.0
    δ_c = binomial(zp.s-1, convert(Int, (zp.s-1)/2))/(2^(zp.s-1))/3 #valore per criterio per α_s 0.15 per 3-mixtures
    #α_start = first(α)
    pars.ψ = 0.99

    results = []
    for β in β
        #### inizializzo vicino alla soluzione a α = 0 ####
        op.m = binomial(zp.s-1,convert(Int, (zp.s-1)/2))/(2^(zp.s-1))#0.5 per 3-mixtures

        for  α in α
            ep.α = α;
            ep.β = β;
    
            println("# NEW ITER: α=$(ep.α)  β=$(ep.β)")
        
            m_old = op.m
            ok, Δ = converge!(op, ep, pars, zp)

            #if α != α_start && op.m < 1e-10
            #    break
            #end

            local δₛ = abs(m_old - op.m)
    
            if δₛ >= δ_c
                open(resfile, "a") do rf
                println(rf, ep.α, " ", 1/ep.β)
                end
    
            push!(results,(ep.α, 1/ep.β))
            break
    
            end
        end
    end

    return results

end