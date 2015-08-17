module GSMDP_mcts

using pattern
using RewardFun
using MCTS_GSMDP


if pattern.nPhases != 1
    error("GSMDP implementation is supposed to work with only 1 phase!")
end

mctsRng = MersenneTwister(myid())

function rollOutPolicy(ST::StateEvent, rngState::AbstractRNG)
    return pattern.g_noaction #our roll-out policy is the silent policy
end

#TODO:
#Since we are using Monte-Carlo, we might be able to use the history
#and actually handle non-exponential time distributions?




function sampleTime(cdf::pattern.cdfTime ,tclip::Float32, p::Float32)
    N = length(cdf.probs)
    
    #adjust to clip time    
    if tclip > 0
        cnt = 0
        for t in cdf.times
            cnt += 1
            if t >= tclip
                break
            end
        end
        p = cdf.probs[cnt] + p * (1f0 - cdf.probs[cnt])
    end

    idx_L = searchsortedlast(cdf.probs, p)
    idx_R = min(idx_L+1, N)

    frac = max((p - cdf.probs[idx_L])/(cdf.probs[idx_R] - cdf.probs[idx_L]), 0.)
    frac = min(frac, 1.)

#    @show idx_L, cdf.probs[idx_L], idx_R, cdf.probs[idx_R],  frac
    t = cdf.times[idx_L] + frac * (cdf.times[idx_R] - cdf.times[idx_L]) 

    return max(t, tclip)
end


function findFirsti(STnow::StateEvent, sojurnTimes::typeof(pattern.sojurnCDFs), rngState::AbstractRNG)
    #TODO: See http://www.ee.ryerson.ca/~courses/ee8103/chap6.pdf
    #for potentiall a better way to select the next state without needing the logs? 
    
    #Find which state will transition
    N = length(Snow)

    tmin = Inf32
    ifirst = 1
    for i in 1:N
       t = sampleTime(pattern.sojurnCDFs[i], STnow[2][i], float32(rand(rngState)))
       if t < tmin
         ifirst = i
         tmin = t::Float32
       end
    end
    
    return (ifirst, tmin)
end


function getNextState!(STnew::StateEvent, STnow::StateEvent, a::typeof(pattern.g_noaction), rngState::AbstractRNG)
        
    #Only transition the one with the earliest event in the race!
    
    
    ifirst = 0;
    t_sojurn = 0f0;
    
    if isNullAct(a)
        (ifirst, t_sojurn) = findFirsti(STnow, pattern.sojurnCDFs, rngState)
    else
        ifirst = a[1];
        t_sojurn = max(pattern.sojurnTime[STnow[1][ifirst]] - STnow[2][ifirst], 0.)
    end
    
    for i in 1:length(Snew)
        if i == ifirst
            STnew[1][ifirst] = randomChoice(STnow[1][ifirst], a[1] == ifirst, a[2], rngState)
            STnew[2][ifirst] = 0f0 #reset event counter since this just transitioned!
        else
            #These did not transition
            STnew[1][i] = STnow[1][i] 
            STnew[2][i] += t_sojurn #increment the sojurn time
        end
    end

    return t_sojurn
end

function getReward(S_T::StateEvent, a::typeof(pattern.g_noaction), pars::MCTS.SPWParams)    
    #assert(pars.β < 0.9f0) #We make the assumption that action cost is small relative to collision cost
    
    #TODO: account for the time!!
    
    S = S_T[1] 
    R = Reward(S, a, pars.β::Float32)
    
    pars.terminate = false;
    #This is a terminal state...
    if( R <=  RewardFun.collisionCost)
        pars.terminate = true
    end
    return R
end

Afun! = pattern.validActions!

assert (typeof(pattern.g_noaction) == MCTS.Action)
assert ((SType, Vector{Float32}) == MCTS.State)

function genMCTSdict(d, ec, n, β, ζ, w, resetDict)
    terminate=false #doesnt matter, getReward will update this at each call
    
    buildTree=false
    pars = MCTS.SPWParams{MCTS.Action}(
                terminate, resetDict, buildTree,    
                d,ec,n,β,ζ, w,
                Afun!,
                rollOutPolicy,
                getNextState!,
                getReward,
                S2LIDX, #unused in GSMP implementation
                mctsRng)
    mcts = MCTS.SPW{MCTS.Action}(pars)
    return mcts
end

###############################
#Default parameters
###############################
d = int16(20*pattern.nPhases)           
ec = abs(RewardFun.collisionCost)*5
n = int32(2000)
β = 0.0f0
ζ = float32(0.5/60.)
w = 5


resetDict = true #reset dictionary every cycle

mcts = genMCTSdict(d, ec, n, β, ζ, w, resetDict)

actWorkspace = Array(extActType, pattern.g_nMaxActs)
actWorkspace[1] = copy(pattern.g_noaction)

function mctsPolicy_gsmdp(S::SType, E::Vector{Float32})
    return MCTS_GSMDP.selectAction!(mcts, actWorkspace, S, E)
end

function loadMCTSPolicy(β::Float32)
    mcts.pars.β = β
    return mctsPolicy_gsmdp
end

export mcts, mctsPolicy_gsmdp, loadMCTSPolicy

# const S = [:LD2, :RB1, :R, :U1]
#  
# function test(S, n)
#      for lo in 1:n 
#          mctsPolicy(S) 
#      end
#  end
#  
# @time test(S,1)
# @time test(S,10)

end
