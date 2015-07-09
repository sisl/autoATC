using MCTS
using pattern

include("CTMDP.jl")

mctsRng = MersenneTwister(1)




function rollOutPolicy(s::SType, rngState)
    return pattern.g_noaction::MCTS.Action #our roll-out policy is the silent policy
end


#TODO: Need to deal with the fact that this is a CTMDP,
#so only one needs to transition at a time!
#Since we are using Monte-Carlo, we might be able to use the history
#and actually handle non-exponential time distributions?
function getNextState(Snow::SType, a::typeof(pattern.g_noaction), rngState)
    Snew = deepcopy(Snow)
    
    #Find which state will transition
    p = rand(rngState, length(Snow))
    t_rand = -Float64[pattern.teaTime[s] for s in Snow] .* log(1-p)
    i = indmin(t_rand)
    
    #Only transition the one with the earliest event in the race!
    Snew[i] = pattern.randomChoice(Snow[i], a[1] == i, a[2]; rngState=rngState)
    
#   This works for an MDP but not a CTMDP !    
#     for i in 1:length(Snow)
#         Snew[i] = pattern.randomChoice(Snow[i], a[1] == i, a[2]; rngState=rngState)
#     end
    return Snew
end

function getReward(S::SType, a::typeof(pattern.g_noaction), pars)
    #TODO: probably can get rid of this
    acomp = pattern.extAct2compAct(a, S)
    
    
    β = pars.β
    assert(β < 0.9f0) #We make the assumption that action cost is small relative to collision cost
    
    
    R = Reward(S, acomp, β)
    
    terminate = false;
    #This is a terminal state...
    if( R <=  collisionCost)
        terminate = true
    end
    return (R, terminate)
end

Afun = pattern.validActions


assert (typeof(pattern.g_noaction) == MCTS.Action)
assert (SType == MCTS.State)




function genMCTSdict(d, ec, n, β)
    pars = MCTS.SPWParams{MCTS.Action}(d,ec,n,β, Afun,rollOutPolicy,getNextState,getReward, mctsRng)
    mcts = MCTS.SPW{MCTS.Action}(pars)
    return mcts
end

d = int16(50)           
ec = abs(collisionCost)
n = int32(1000)
β = 0.0f0
mcts = genMCTSdict(d, ec, n, β)


# 
# gridworld = GenerativeModel(getInitialState,getNextState,getReward)
# 
# simulate(gridworld,mcts,policy,nSteps,rng)