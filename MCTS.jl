module MCTS
# This module implements Monte Carlo tree search (MCTS), an online MDP solver.
#It is adapted by <zouhair> from code available @ https://bitbucket.org/sisl/mdp.jl/


# In MCTS. the poosible evolutions of the system are represented as a tree
# and an approximation to this tree is built iteratively and used to inform the choice of action.
#Documentation for MCTS can be found in Decision Making Under Uncertainty, MIT Press 2015,
#and A survey of Monte Carlo tree search methods in Computational Intelligence and AI in games, 2012.

export SPWParams, SPW, selectAction!

typealias Depth Int16
typealias Reward Float32


#TODO: This forces a certain type for state and action ...
typealias State Array{Symbol,1}
typealias Action (Int8, Symbol)

type SPWParams{T<:Action}
    d::Depth                    # search depth
    ec::Float32                 # exploration constant- governs trade-off between exploration and exploitation in MCTS
    n::Int32                    # number of iterations
    β::Float32
    
    A::Function                 # set of allowable actions 
    rolloutPolicy::Function     # returns action for rollout policy
    getNextState::Function      # takes state and action as arguments and returns next state from generative model
    getReward::Function         # takes state and action as arguments and returns reward
    
    rng::AbstractRNG            # random number generator

end


#statistics for a given node
type StateStat
    n::Array{Int32,1}  #Number of times the state/action pairs were visitied
    q::Array{Reward,1} #average reward for that state
    
    StateStat(nA) = new(zeros(Int32,nA),zeros(Reward,nA))
end

type SPW{T}
    stats::Dict{State,StateStat} #statistics
    pars::SPWParams{T}          #parameters
    
    #constructor (takes parameters, initializes statistics to be empty)
    SPW{T<:Action}(p::SPWParams{T}) = new(Dict{State,StateStat}(),p)
end

#This is the entry function, it needs an initial (root) state, and
#the parameters for the SPW structure which contains the parameters
#and the statistics for the tree
# This function calls simulate and chooses the approximate best action from the reward approximations 
function selectAction!(spw::SPW,s0::State)
    
    #TODO: try to call A as little as possible?
    acts = spw.pars.A(s0) #get the allowable actions
    
    
    #This is to avoid the first call to simulate wasting a rollout
    if !haskey(spw.stats, s0)
        spw.stats[s0] = StateStat(length(acts))
    end
    
    for i = 1:spw.pars.n 
        simulate(spw,s0,spw.pars.d)
    end
    
    #Choose action with highest apporoximate q-value 
    return acts[indmax(spw.stats[s0].q)]::Action 
end


#TODO: handle terminal states!
function simulate(spw::SPW,s::State,d::Depth)
    # This function returns the reward for one iteration of MCTS
    
    #We have reached the bottom, bubble up.
    if d == 0
        return 0.0f0::Reward
    end
    
    #Determine actions available for this state
    acts = spw.pars.A(s)
    nActs = length(acts)
    
    #If this state has no statistics yet (i.e. first visit)
    #perform a roll-out
    
    
    if !haskey(spw.stats,s)
        #Add a new node with the state statistics
        spw.stats[s] = StateStat(nActs)
        #perform one rollout and return
        return rollout(spw,s,d)::Reward
    else # choose an action using UCT criterion
        
        
        # choose action with highest UCT score
        # which trades-off exploration / explotation
        cS = spw.stats[s]
        
        N = sum(cS.n)
        #Check N > 0 to avoid log(0)/0 = -Inf, can't take sqrt!
        i = 1 #if N == 0, use first action
        if( N > 0)
            logN =  log(N) + eps(0.) # +eps is to avoid log(1)/0 = NaN
            i  = indmax(cS.q + spw.pars.ec * sqrt( logN ./ cS.n))
        end
        a  = acts[i]
        
        #Randomly select the next state based on the action used
        sp = spw.pars.getNextState(s,a,spw.pars.rng)
        #Estimate the reward
        (q, terminate) = spw.pars.getReward(s,a, spw.pars) 
        if !terminate
            q += simulate(spw,sp,int16(d-1))
        end
        
        #Update the statistics
        cS.n[i] += one(Int32)
 
        #This could maybe take into account the transition probablities?
        #i.e. given that we know P(s' | s, a), use importance weighing?
        #Doh, apparently people thought about this already. See:
        #http://www.ai.rug.nl/~mwiering/Tom_van_der_Kleij_Thesis.pdf
        #Equ 3.2
        cS.q[i] += (q-cS.q[i])/cS.n[i]
        return q::Reward
    end
end

function rollout(spw::SPW,s::State,d::Depth)
    # Runs a rollout simulation using the default policy
    if d == 0
        return 0.0f0::Reward
    else 
        a  = spw.pars.rolloutPolicy(s,spw.pars.rng)
        sp = spw.pars.getNextState(s,a,spw.pars.rng)
        #This runs a roll-out simulation, note the lack of discount factor...
        (R, terminate) = spw.pars.getReward(s,a, spw.pars)
        if !terminate
            R +=  rollout(spw,sp,int16(d-1))::Reward
        end
        return R::Reward
    end 
end

end # module