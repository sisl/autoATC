using MCTS
using pattern

include("CTMDP.jl")

mctsRng = MersenneTwister(1)




function rollOutPolicy(s::SType, rngState::AbstractRNG)
    return pattern.g_noaction::MCTS.Action #our roll-out policy is the silent policy
end


#TODO: Need to deal with the fact that this is a CTMDP,
#so only one needs to transition at a time!
#Since we are using Monte-Carlo, we might be able to use the history
#and actually handle non-exponential time distributions?
function findFirsti(Snow::SType, trantimes::typeof(pattern.teaTime), rngState::AbstractRNG)
    #Find which state will transition
    N = length(Snow)

    tmin = Inf
    ifirst = 1
    for i in 1:N
       t = -trantimes[Snow[i]] * log(1-rand(rngState))
       if t < tmin
         ifirst = i
         tmin = t
       end
    end
    
    return ifirst
end



function getNextState!(Snew::SType, Snow::SType, a::typeof(pattern.g_noaction), rngState::AbstractRNG)
    #copy!(Snew, Snow)
    
    #Only transition the one with the earliest event in the race!
    ifirst = findFirsti(Snow, pattern.teaTime, rngState)
    #Snew[ifirst] = randomChoice(Snow[ifirst], a[1] == ifirst, a[2], rngState)
    
    for i in 1:length(Snew)
        if i == ifirst
            Snew[ifirst] = randomChoice(Snow[ifirst], a[1] == ifirst, a[2], rngState)
        else
            Snew[i] = Snow[i] 
        end
    end
#     
    return nothing
end

function getReward(S::SType, a::typeof(pattern.g_noaction), pars::MCTS.SPWParams)    
    #assert(pars.β < 0.9f0) #We make the assumption that action cost is small relative to collision cost
    
    R = Reward(S, a, pars.β::Float32)
    
    pars.terminate = false;
    #This is a terminal state...
    if( R <=  collisionCost)
        pars.terminate = true
    end
    return R
end

Afun! = pattern.validActions!


assert (typeof(pattern.g_noaction) == MCTS.Action)
assert (SType == MCTS.State)




function genMCTSdict(d, ec, n, β, resetDict )
    terminate=false#doesnt matter, getReward will update this at each call
    
    pars = MCTS.SPWParams{MCTS.Action}(terminate, resetDict, d,ec,n,β, Afun!,rollOutPolicy,getNextState!,getReward, S2LIDX, mctsRng)
    mcts = MCTS.SPW{MCTS.Action}(pars)
    return mcts
end

d = int16(50)           
ec = abs(collisionCost)
n = int32(1000)
β = 0.0f0
resetDict = true #reset dictionary every cycle

mcts = genMCTSdict(d, ec, n, β, resetDict)

actWorkspace = Array(extActType, g_nCompActs)
actWorkspace[1] = copy(pattern.g_noaction)

mctsPolicy = S -> MCTS.selectAction!(mcts, actWorkspace, S)




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
