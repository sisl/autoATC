module GSMDP_mcts

using pattern
using RewardFun
using MCTS_GSMDP

typealias StateEvent MCTS_GSMDP.StateEvent

if pattern.nPhases != 1
    error("GSMDP implementation is supposed to work with only 1 phase!")
end

mctsRng = MersenneTwister(myid())

function rollOutPolicy(ST::StateEvent, rngState::AbstractRNG)
    return pattern.g_noaction #our roll-out policy is the silent policy
end



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

    frac = max((p - cdf.probs[idx_L])/(cdf.probs[idx_R] - cdf.probs[idx_L]), 0)
    frac = min(frac, 1.)

#    @show idx_L, cdf.probs[idx_L], idx_R, cdf.probs[idx_R],  frac
    t = cdf.times[idx_L] + frac * (cdf.times[idx_R] - cdf.times[idx_L]) 

    return float32(max(t, tclip))
end


function findFirsti(STnow::StateEvent, sojurnCDF::typeof(pattern.sojurnCDF), rngState::AbstractRNG)
    #TODO: See http://www.ee.ryerson.ca/~courses/ee8103/chap6.pdf
    #for potentiall a better way to select the next state without needing the logs? 
    
    #Find which state will transition
    N = length(STnow[1])

    tmin = Inf32
    ifirst = 1
    for i in 1:N
       t = sampleTime(sojurnCDF[STnow[1][i]], STnow[2][i], float32(rand(rngState))) - STnow[2][i]
       if t < tmin
         ifirst = i
         tmin = t::Float32
       end
    end
    
    #round to 5 seconds ...
    tmin = round(tmin/5)*5
    
    return (ifirst, tmin)
end


function getNextState!(STnew::StateEvent, STnow::StateEvent, anew::typeof(pattern.g_noaction), rngState::AbstractRNG)
    #Only transition the one with the earliest event in the race!
    
    ifirst = 0;
    t_sojurn = 0f0;
    
    activeAct = pattern.compAct2extAct(STnow[3], STnow[1])
    if pattern.isNullAct(activeAct)
        activeAct = anew
    elseif !pattern.isNullAct(anew)
        error("Cannot issue action if the previous action has not been acted upon!")
    end


    (ifirst, t_sojurn) = findFirsti(STnow, pattern.sojurnCDF, rngState)
    
    for i in 1:length(STnew[1])
        if i == ifirst
            STnew[1][ifirst] = randomChoice(STnow[1][ifirst], activeAct[1] == ifirst, activeAct[2], rngState)
            STnew[2][ifirst] = 0f0 #reset event counter since this just transitioned!
        else
            #These did not transition
            STnew[1][i] = STnow[1][i] 
            STnew[2][i] = STnow[2][i] + t_sojurn #increment the sojurn time
        end
    end


    #Assume action was acted upon
    STnew[3][1] = 0
    STnew[3][2] = 0
    
    if !pattern.isNullAct(activeAct)
        if ifirst != activeAct[1]
            compact = pattern.extAct2compAct(activeAct, STnew[1])
            STnew[3][1] = compact[1]
            STnew[3][2] = compact[2]
        end
    end
    return t_sojurn
end


cnt = 10
function getReward(S_T::StateEvent, a::typeof(pattern.g_noaction), pars::MCTS_GSMDP.SPWParams)    
    #assert(pars.β < 0.9f0) #We make the assumption that action cost is small relative to collision cost
    
    S = S_T[1] 
    R = Reward(S, a, pars.β::Float32, timeHorizon = 5f0, E = S_T[2])
    
    
    pars.terminate = false;
    #This is a terminal state...
    if( R <=  RewardFun.collisionCost)
        pars.terminate = true
#         global cnt
#         if cnt > 0
#             println("Terminated @ $S_T")
#             cnt -= 1
#         end
    end
    return R
end

function Afun!(acts::Vector{MCTS_GSMDP.Action}, S_T::StateEvent)

    nActs = pattern.validActions!(acts, S_T[1])
    #restrict to null act if action has not yet been executed!
    if !pattern.isNullAct(S_T[3])
        nActs = 1
    end
    
    return nActs
     
end

assert (typeof(pattern.g_noaction) == MCTS_GSMDP.Action)
assert ((SType, Vector{Float32}, typeof(pattern.g_nullAct)) == MCTS_GSMDP.StateEvent)
assert (SType == MCTS_GSMDP.State)

function genMCTSdict(d, ec, n, β, ζ, w, resetDict)
    terminate=false #doesnt matter, getReward will update this at each call
    
    buildTree=false
    pars = MCTS_GSMDP.SPWParams{MCTS_GSMDP.Action}(
                terminate, resetDict, buildTree,    
                d,ec,n,β,ζ, w,
                Afun!,
                rollOutPolicy,
                getNextState!,
                getReward,
                S2LIDX, #unused in GSMP implementation
                mctsRng)
    mcts = MCTS_GSMDP.SPW{MCTS_GSMDP.Action}(pars)
    return mcts
end

###############################
#Default parameters
###############################
d = int16(200)           
ec = abs(RewardFun.collisionCost)*5
n = int32(2000)
β = 0.0f0
ζ = float32(0.5/60.)
w = 100


resetDict = true #reset dictionary every cycle

mcts = genMCTSdict(d, ec, n, β, ζ, w, resetDict)

actWorkspace = Array(extActType, pattern.g_nMaxActs)
actWorkspace[1] = copy(pattern.g_noaction)

function mctsPolicy_gsmdp(S::SType, E::Vector{Float32})
    return MCTS_GSMDP.selectAction!(mcts, actWorkspace, (S,E,pattern.g_nullAct))
end

function loadMCTSPolicy_gsmdp(β::Float32)
    mcts.pars.β = β
    return mctsPolicy_gsmdp
end

export mcts, mctsPolicy_gsmdp, loadMCTSPolicy_gsmdp

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
