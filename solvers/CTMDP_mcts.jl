module CTMDP_mcts

using pattern
using RewardFun
using MCTS

using CTMDP_kronindexing

mctsRng = MersenneTwister(myid())

function rollOutPolicy(s::SType, rngState::AbstractRNG)
    return pattern.g_noaction::MCTS.Action #our roll-out policy is the silent policy
end

#TODO:
#Since we are using Monte-Carlo, we might be able to use the history
#and actually handle non-exponential time distributions?
function findFirsti(Snow::SType, sojurnTimes::typeof(pattern.sojurnTime), rngState::AbstractRNG)
    #TODO: See http://www.ee.ryerson.ca/~courses/ee8103/chap6.pdf
    #for potentiall a better way to select the next state without needing the logs? 
    
    #Find which state will transition
    N = length(Snow)

    tmin = Inf32
    ifirst = 1
    for i in 1:N
       t = -sojurnTimes[Snow[i]] * log(1-float32(rand(rngState)))
       if t < tmin
         ifirst = i
         tmin = t::Float32
       end
    end
    
    return (ifirst, tmin)
end


function getNextState!(Snew::SType, Snow::SType, a::typeof(pattern.g_noaction), rngState::AbstractRNG)
    #Only transition the one with the earliest event in the race!
    
    
    ifirst = 0;
    t_sojurn = 0f0;
    
    if isNullAct(a)
        (ifirst, t_sojurn) = findFirsti(Snow, pattern.sojurnTime, rngState)
    else
        ifirst = a[1];
        t_sojurn = pattern.sojurnTime[Snow[ifirst]]
    end
    
    for i in 1:length(Snew)
        if i == ifirst
            Snew[ifirst] = randomChoice(Snow[ifirst], a[1] == ifirst, a[2], rngState)
        else
            Snew[i] = Snow[i] 
        end
    end
#     
    return t_sojurn
end

function getReward(S::SType, a::typeof(pattern.g_noaction), pars::MCTS.SPWParams)    
    #assert(pars.β < 0.9f0) #We make the assumption that action cost is small relative to collision cost
    
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
assert (SType == MCTS.State)

function genMCTSdict(d, ec, n, β, ζ, resetDict)
    terminate=false#doesnt matter, getReward will update this at each call
    
    buildTree=false
    pars = MCTS.SPWParams{MCTS.Action}(
                terminate, resetDict, buildTree,    
                d,ec,n,β,ζ, 
                Afun!,
                rollOutPolicy,
                getNextState!,
                getReward,
                S2LIDX,
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

resetDict = true #reset dictionary every cycle

mcts = genMCTSdict(d, ec, n, β, ζ, resetDict)

actWorkspace = Array(extActType, pattern.g_nMaxActs)
actWorkspace[1] = copy(pattern.g_noaction)

function mctsPolicy(S::SType, E::Vector{Float32})  
    return MCTS.selectAction!(mcts, actWorkspace, S)
end

function loadMCTSPolicy(β::Float32)
    mcts.pars.β = β
    return mctsPolicy
end

export mcts, mctsPolicy, loadMCTSPolicy

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
