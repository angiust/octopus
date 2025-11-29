module P

using QuadGK
using ForwardDiff
using Optim

include("../common/common.jl")


###### INTEGRATION  ######
const ∞ = 10.0
const dx = 0.01

const interval = map(x->sign(x)*abs(x)^2, -1:dx:1) .* ∞

∫D(f, int=interval) = quadgk(z->begin
        r = G(z) .* f(z)
        isfinite(r) ? r : 0.0
    end, int..., atol=1e-7, maxevals=10^7)[1]


############### PARAMS ################

@with_kw mutable struct OrderParams
	q::Float64 = 0.5
    qh::Float64 = 0.2
    m::Float64 = 0.3
end

collect(op::OrderParams) = [getfield(op, f) for f in fieldnames(op)]

@with_kw mutable struct ExtParams
    α::Float64 = 0.1  # constrained density
    β::Float64 = 1.0  # inverse temperature
    d::Float64 = 0.0 # dilution
end


@with_kw mutable struct Params
    ϵ::Float64 = 1e-5       # stop criterium
    ψ::Float64 = 0.9         # damping
    maxiters::Int = 10000
    verb::Int = 0
end

mutable struct ThermFunc
    f::Float64 		# free energy
end


Base.show(io::IO, op::OrderParams) = shortshow(io, op)
Base.show(io::IO, ep::ExtParams) = shortshow(io, ep)
Base.show(io::IO, tf::ThermFunc) = shortshow(io, tf)

###################################################################################

#### INTERACTION TERM ####
function Gi(q, qh, m, α, β)
    α / 2 + 0.5 * m^2 + 0.5 * α * β * qh * (1 - q)
end

Gi(op::OrderParams, ep::ExtParams) = Gi(op.q, op.qh, op.m, ep.α, ep.β)

function Gs(q, α, β, d)
    term = 1 - β * (1-d) * (1-q)
    0.5 * α * (log(term) - β * (1-d) * q / term ) / β
end

Gs(op::OrderParams, ep::ExtParams) = Gs(op.q, ep.α, ep.β, ep.d)


#### ENERGETIC TERM ####
function Ge(qh, m, α, β, d)
    active_term = ∫D(z -> begin
        term = √(α * qh) * z + m
        log(2 * cosh(β * term))
    end)

    inactive_term = ∫D(z -> begin
        term = √(α * qh) * z
        log(2 * cosh(β * term))
    end)

    return - ((1-d) * active_term + d * inactive_term) / β
end

Ge(op::OrderParams, ep::ExtParams) = Ge(op.qh, op.m, ep.α, ep.β, ep.d)

#### SADDLE POINT EQUATIONS SINGLE PATTERN RETRIEVAL ####
function update_qh(q, β, d)
    (1 - d)^2 * q / (1 - β * (1 - d) * (1 - q))^2
end

update_qh(op::OrderParams, ep::ExtParams) = update_qh(op.q, ep.β, ep.d)


function update_q(qh, m, α, β, d)
    active_term = ∫D(z -> begin
        term = √(α * qh) * z + m
        tanh(β * term)^2
    end)
    inactive_term = ∫D(z -> begin
        term = √(α * qh) * z
        tanh(β * term)^2
    end)

    return (1-d) * active_term + d * inactive_term
end

update_q(op::OrderParams, ep::ExtParams) = update_q(op.qh, op.m, ep.α, ep.β, ep.d)

function update_m(qh, m, α, β, d)
    (1 - d) * ∫D(z -> begin
        term = √(α * qh) * z + m
        tanh(β * term)
    end)
end

update_m(op::OrderParams, ep::ExtParams) = update_m(op.qh, op.m, ep.α, ep.β, ep.d)



############ Thermodynamic functions ############

function free_energy(op::OrderParams, ep::ExtParams)
    Gi(op, ep) + Gs(op, ep) + Ge(op, ep)
end

## Thermodynamic functions
function all_therm_func(op::OrderParams, ep::ExtParams)
    f = free_energy(op, ep)

    return ThermFunc(f)
end

#################  SADDLE POINT  ##################
# Right-hand-side

fq(op, ep) = update_q(op, ep)			   	   # q = fq (der: qh)
fqh(op, ep) = update_qh(op, ep)		   		   # qh = fqh (der: q)
fm(op, ep) = update_m(op, ep)			   	   # m = fm (der: m)


function converge!(op::OrderParams, ep::ExtParams, pars::Params)
    @extract pars: maxiters verb ϵ ψ
    Δ = Inf
    ok = false   	

    for it = 1:maxiters
        Δ = 0.0
        ok = true
        verb > 1 && println("########## it=$it ##########")

        ########################################################################
        local op1 = deepcopy(op)

        @update  op.q    fq       Δ ψ verb  op1 ep   #update q
        @update  op.qh   fqh       Δ ψ verb  op1 ep  #update qh
        @update  op.m    fm       Δ ψ verb  op1 ep   #update m
        
        ########################################################################

        verb > 1 && println(" Δ=$Δ\n")
        #(println(op); println(ep))

        @assert isfinite(Δ)
        ok &= Δ < ϵ
        ok && break         # if ok==true, exit
    end
    ok, Δ
end


function converge(;
    q = 0.5, qh = 0.2, m = 0.3,
    α = 0.1, β = 1.0, d = 0.0,
    ϵ = 1e-5, maxiters = 10000, verb = 3, ψ = 0.,
    kws...)

    op = OrderParams(q, qh, m)
    ep = ExtParams(α, β, d)
    pars = Params(ϵ, ψ, maxiters, verb)
	
    converge!(op, ep, pars)
    tf = all_therm_func(op, ep)
    println(tf)
    return op, ep, pars, tf
end



"""
    readparams(file::String, line::Int=-1)

Read order and external params from results file.
Zero or negative line numbers are counted
from the end of the file.
"""
function readparams(file::String, line::Int=0)
    data = readdlm(file, String)
    l = line > 0 ? line : length(data[:,2]) + line
    v = map(x-> begin
                    try
                        parse(Float64, x)
                    catch
                        Symbol(x)
                    end
                end, data[l,:])
    #return v
    i0 = length(fieldnames(ExtParams))
    i1 = i0 + 1  + length(fieldnames(ThermFunc))
    iend = i1 - 1 + length(fieldnames(OrderParams))

    return ExtParams(v[1:i0]...), OrderParams(v[i1:iend]...)
end


function span(;
        q = 0.5, qh = 0.2, m = 0.3,
        α = 0.1, β = 1.0, d = 0.0,
        ϵ = 1e-5, maxiters = 10000, verb = 3, ψ = 0.,
        kws...)

    op = OrderParams(q, qh, m)
    ep = ExtParams(first(α), first(β), d) # attention: i think d is fixed, is it right?
    pars = Params(ϵ, ψ, maxiters, verb)

    return span!(op, ep, pars; α=α, β=β, kws...) # should i pass d explicitly too? for now i keep it fixed so it should be ok
end

function span!(op::OrderParams, ep::ExtParams, pars::Params;
        α=0.2, β=1.0, 
        resfile = "hopfield.txt")

    !isfile(resfile) && open(resfile, "w") do f
        allheadersshow(f, ExtParams, ThermFunc, OrderParams)
    end

    results = []

    for β in β, α in α
        ep.α = α;
        ep.β = β;

        println("# NEW ITER: α=$(ep.α)  β=$(ep.β)")
        
	    ok, Δ = converge!(op, ep, pars)

        tf = all_therm_func(op, ep)

        push!(results, (ok, deepcopy(ep), deepcopy(tf), deepcopy(op)))

        ok && open(resfile, "a") do rf
            println(rf, plainshow(ep), " ", plainshow(tf), " ", plainshow(op))
        end

        !ok && break # a α = 0.138 smette di convergere e quindi fa il break
        pars.verb > 0 && print(ep, "\n", tf,"\n")
    end
    return results
end


end #module