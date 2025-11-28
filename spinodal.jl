################################################
################################################
#
#
#   The spinodal lines are given by the stability 
#   matrix's lowest eigenvalue, the condition is 
#   λ=0 and it can be solved using newton method
#
################################################
################################################


using Parameters
using ExtractMacro


include("../common/common.jl")
include("./hopfield.jl")
include("mixtures.jl")



#################   LAMBDA STABILITY MIXTURES   ##################
function get_q_bar(m, qh, β, α, d, s)
    # when s = 1 is different because we don't have the distribution ζ_{s-1} and the average respects to it
    if s == 1
        # A = √(α*qh)*z. avg_xi_mu = 0.5 * [tanh²(A+m)+tanh²(A-m)]
        avg_val = P.∫D(z -> begin
            A = √(α * qh) * z
            0.5 * (tanh(β * (A + m))^2 + tanh(β * (A - m))^2)
        end)
        return (1 - d) * avg_val
    end    
    # if s > 1 we need the distribution ζ_{s-1}
    zp_n_minus_1 = calc_zeta_prob(s-1, d)
    zeta_nm1 = zp_n_minus_1.ζ
    prob_nm1 = zp_n_minus_1.prob

    # two (three with z) level average
    avg = P.∫D(z -> begin
        result = 0.0
        # average on ζ_{n-1}
        for i in eachindex(zeta_nm1)
            A = √(α * qh) * z + m * zeta_nm1[i]
            # inner average on ξ^μ
            avg_xi_mu = 0.5 * (tanh(β * (A + m))^2 + tanh(β * (A - m))^2)
            result += prob_nm1[i] * avg_xi_mu
        end
        return result
    end)

    (1-d) * avg
end
get_q_bar(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = get_q_bar(op.m, op.qh, ep.β, ep.α, ep.d, zp.s)

function get_Q(m, q, qh, β, α, d, s, ζ, prob)
    # s = ζ[end]
    @assert s > 1 "Q is defined only for s > 1"

    term1 = P.∫D(z -> begin
        result = 0.0
        for i in eachindex(ζ) #### valore medio sulla distribuzione di ζ
            term = √(α * qh) * z + m * ζ[i]
            result += ζ[i]^2 * tanh(β * term)^2 * prob[i]
        end
        result # now i don't need to multiply by 2 because the sum is over all the possible values of ζ
    end)

    q_bar = get_q_bar(m, qh, β, α, d, s)

    term1 / (s * (s - 1)) - q_bar / (s-1)
end
get_Q(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = get_Q(op.m, op.q, op.qh, ep.β, ep.α, ep.d, zp.s, zp.ζ, zp.prob)

function get_λ3(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb, Q::AbstractFloat, q_bar::AbstractFloat)
    1.  - ep.β * (1 - ep.d - q_bar) - ep.β * Q
end
get_λ3(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = get_λ3(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb, get_Q(op, ep, zp), get_q_bar(op, ep, zp))

function get_λ1(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb, Q::AbstractFloat, q_bar::AbstractFloat)
    1.  - ep.β * (1 - ep.d - q_bar) + (zp.s-1)ep.β * Q
end
get_λ1(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = get_λ1(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb, get_Q(op, ep, zp), get_q_bar(op, ep, zp))

function get_λ2(op::P.OrderParams, ep::P.ExtParams)
    1. - ep.β * (1 - ep.d) * (1 - op.q)
end

# for now i ignore λrsb, which i do not kwown how it is with dilution
function get_λrsb(m,q,qh,α,β,ζ,prob)
    s = ζ[end]        
    P.∫D(z -> begin
        result = 0.0
        for i in eachindex(ζ) #### valore medio sulla distribuzione di ζ
            term = √(α * qh) * z + m * ζ[i] 
            result += sech(β*term)^4 * prob[i]
        end
        (1-β*(1-q))^2-α*β*β*2*result
    end)
end
get_λrsb(op::P.OrderParams, ep::P.ExtParams, zp::ZetaProb) = get_λrsb(op.m,op.q,op.qh,ep.α,ep.β,zp.ζ,zp.prob)

function find_spinodal(n_points::Int, s::Int = 3, resfile = "prova.txt"; atol = 1e-6)
    !isfile(resfile) && open(resfile, "w") do rf
        println(rf,"# spinodal line for $s-mixtures")
        end

    zp = calc_zeta_prob(s)
    
    T = 0.05:(0.455-0.05)/n_points:0.455
    β = map(x -> 1. /x, reverse(T))
    
    δα = 1e-5

    op = P.OrderParams(0.5, 1., 0.3)
    ep = P.ExtParams(0.0, first(β))
    pars = P.Params(ψ = 0.5)

    results = []
    λ = Inf

    for β in β
        # a β fissata facciamo un while su α fino a che non troviamo λ < atol, opppure lo superiamo,
        # in quel caso torniamo indietro e dimezziamo il passo in α
        it = 0
        ep.β = β
        println("# β = $β")
        while true
            
            it += 1
            converge!(op, ep, pars, zp)
            λ = get_λ(op, ep, zp) #calcola l'autovalore
            if λ < 0
                # facciamo un passo indietro e dimezziamo il passo in α fino a che λ > 0
                while λ < 0
                    ep.α -= δα  
                    converge!(op, ep, pars, zp)
                    δα = δα /2.
                    ep.α += δα
                    converge!(op, ep, pars, zp)

                end
            end

            if λ >= 0 && λ < atol
                #se troviamo α tale che λ è sotto la soglia di tolleranza, abbiamo finito

                break;
            end
            

            println("it: $it  m: $(op.m)  q: $(op.q)  r: $(op.qh)   λ: $λ   α: $(ep.α)")
            ep.α +=δα
        end
        println("λ: $λ α: $(ep.α)   β: $(ep.β)")
        push!(results,(λ, ep.α, ep.β))
        

        open(resfile,"a") do rf
            println(rf,"λ: $λ α: $(ep.α)   β: $(ep.β)")
        end

        ep.α -=  δα
    end

    return results

end

function mix_span(αn::Array{Float64}, βn::Array{Float64}, s::Int = 3, d::Float64 = 0.0, resfile = "./data_3.txt")
    resfile = "result/mix_s=$(s)_d=$(d).txt"

    !isfile(resfile) && open(resfile, "w") do rf
        println(rf,"# $s-mixtures   1:α 2:T 3:m 4:q 5:r 6:Q 7:q_bar 8:λ3 9:λ1 10:λ2")
        end

        op = P.OrderParams(0.3, 1, 0.5)
        ep = P.ExtParams()
        ep.d = d
        pars = P.Params(ψ = 0.8)

        zp = calc_zeta_prob(s, d)
        for β ∈ βn
            ep.β = β
            println("#iteration β = $β")
            op.m = 0.5
            for α in αn
                ep.α = α
                converge!(op, ep, pars, zp)
                Q = get_Q(op, ep, zp)
                q_bar = get_q_bar(op, ep, zp)
                λ3 = get_λ3(op, ep, zp, Q, q_bar)
                λ1 = get_λ1(op, ep, zp, Q, q_bar)
                λ2 = get_λ2(op, ep)
                open(resfile, "a") do rf
                    println(rf,"$(ep.α) $(1/ep.β) $(op.m) $(op.q) $(op.qh) $Q $q_bar $(λ3) $(λ1) $(λ2)")
                    end
            end
        end

    end



    #= function generate_patterns(N::Int)
        if N == 1
            return [[1], [-1]]
        else
            sub_vectors = generate_patterns(N-1)
            vectors = []
            for v in sub_vectors
                push!(vectors, [1; v])
                push!(vectors, [-1; v])
            end
            return vectors
        end
    end =#