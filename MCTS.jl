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
    terminate::Bool             #Flag to be populated by getReward to indicate we have reached terminal state
    resetDict::Bool             #Whether to reset dictionary
    
    d::Depth                    # search depth
    ec::Float32                 # exploration constant- governs trade-off between exploration and exploitation in MCTS
    n::Int32                    # number of iterations
    β::Float32                  # Alert/Collision ration (should be inside of problem defintion...)
    γ::Float32                  # Discount factor
    
    Afun!::Function              # set of allowable actions 
    rolloutPolicy::Function     # returns action for rollout policy
    getNextState!::Function     # takes state and action as arguments and populates next state from generative model
    getReward::Function         # takes state and action as arguments and returns reward
    hashState::Function
    rng::AbstractRNG            # random number generator

end


#statistics for a given node
type StateStat
    n::Array{Float32,1}  #Number of times the state/action pairs were visitied
    q::Array{Reward,1} #average reward for that state
    
    StateStat(nA) = new(zeros(Float32,nA),zeros(Reward,nA))
end



#Using our own hash to avoid weirdness with dictionaries...
typealias StateKey Int64
type SPW{T}
    stats::Dict{StateKey,StateStat} #statistics
    pars::SPWParams{T}          #parameters
    
    #constructor (takes parameters, initializes statistics to be empty)
    #Note that we are 
    SPW{T<:Action}(p::SPWParams{T}) = new(Dict{StateKey,StateStat}(),p)
end

#This is the entry function, it needs an initial (root) state, and
#the parameters for the SPW structure which contains the parameters
#and the statistics for the tree
# This function calls simulate and chooses the approximate best action from the reward approximations 
function selectAction!(spw::SPW, acts::Vector{Action}, s0::State)
       
    nActs = spw.pars.Afun!(acts, s0)
    
    #If there's only one action to be taken, no point in wasting time...
    if(nActs == 1)
        return acts[1]
    end

    #Reset the dictionary if told to
    if(spw.pars.resetDict)
        #TODO: compare the performance to using empty!(spw.stats) ?
        #spw.stats = Dict{StateKey,StateStat}()
        for v in values(spw.stats)
            fill!(v.n, 0.0f0)
            fill!(v.q, 0.0f0)
        end
    end    
    
    #This is to avoid the first call to simulate wasting a rollout
    s0Hash = spw.pars.hashState(s0)::StateKey
    if !haskey(spw.stats, s0Hash)
        spw.stats[s0Hash] = StateStat(nActs)
    end



    #TODO: allocate these outside and pass them afterward
    #we need this to be able to restore the state!    
    s0_orig = deepcopy(s0)
    #This is a workspace variable to avoid having to re-allocate
    sp0 = deepcopy(s0)
    
    for i = 1:spw.pars.n 
        simulate!(spw,acts, s0,sp0,spw.pars.d)
        copy!(s0, s0_orig) #make sure we reset s0 everytime!
    end

    #get the list of allowable actions for this state
    #We need to do this since simulate! will work on acts
    
    #TODO: Find a way to avoid calling this function a 2nd time!
    
    nActs = spw.pars.Afun!(acts, s0)     
    #Choose action with highest apporoximate q-value 
    return acts[indmax(spw.stats[s0Hash].q)]::Action 
end



#Sp is passed-in to avoid to have to allocate memory
#Note that this function will mangle s
function simulate!(spw::SPW, acts::Vector{Action}, s::State, sp::State, d::Depth)
    # This function returns the reward for one iteration of MCTS
    
    #We have reached the bottom, bubble up.
    if d == 0
        return 0.0f0::Reward
    end
    
#      tabs = string([" " for i in 1:(2*d)]...)

    
    #println(tabs, "(",d,")", "Entry s=",s, " -> sp =",sp)
    
    #TODO: Make the A(s) function better!
    #Determine actions available for this state
    nActs = spw.pars.Afun!(acts, s)
    #nActs = length(acts)::Int64
    
    #If this state has no statistics yet (i.e. first visit)
    #perform a roll-out
    
    
    skey = spw.pars.hashState(s)::StateKey
    if !haskey(spw.stats,skey)
        #Add a new node with the state statistics
        spw.stats[skey] = StateStat(nActs)
        #perform one rollout and return
        return rollout!(spw,s,sp,d)::Reward
    else # choose an action using UCT criterion
        
        
        # choose action with highest UCT score
        # which trades-off exploration / explotation
        cS = spw.stats[skey]
        
        N = sum(cS.n)::Float32
        #Check N > 0 to avoid log(0)/0 = -Inf, can't take sqrt!
        i = 1 #if N == 0, use first action
        if( N > 0)
            ec_sqrt_logN =  (spw.pars.ec * sqrt(log(N) + eps(0.0f0)))::Float32 # +eps is to avoid log(1)/0 = NaN
            
            #i  = indmax(cS.q + ec_sqrt_logN ./ sqrt(cS.n))
            maxQec = -Inf32
            q_ec = 0.0f0
            for j = 1:nActs
                q_ec = (cS.q[j] +  ec_sqrt_logN / sqrt(cS.n[j]))::Reward
                if q_ec > maxQec
                    i = j
                    maxQec = q_ec
                end
            end
        end
        a  = acts[i]
        
        #Randomly select the next state based on the action used
        spw.pars.getNextState!(sp, s, a,spw.pars.rng)
        #println(tabs, "(",d,")", "Next  s=",s, " -> sp =",sp)
        #Estimate the reward
        q = spw.pars.getReward(s,a, spw.pars)::Reward 
        if !spw.pars.terminate
            #Note that a call to simulate! will change s and sp, but we 
            #don't care at this point since we no longer need their values!
            q += spw.pars.γ*simulate!(spw,acts, sp,s,int16(d-1)) #Note that s is now just a temp storage variable
        end
        #println(tabs, "(",d,")", "Exit  s=",s, " -> sp =",sp)

        #Update the statistics
        cS.n[i] += 1.0f0
 
        #This could maybe take into account the transition probablities?
        #i.e. given that we know P(s' | s, a), use importance weighing?
        #Doh, apparently people thought about this already. See:
        #http://www.ai.rug.nl/~mwiering/Tom_van_der_Kleij_Thesis.pdf
        #Equ 3.2
        cS.q[i] += (q-cS.q[i])/cS.n[i]
        return q::Reward
    end
end

function rollout!(spw::SPW,s::State,sp::State,d::Depth)
    # Runs a rollout simulation using the default policy
#      tabs = string([" " for i in 1:(2*d)]...)

    R = 0.0f0::Reward
    if d > 0
        #println(tabs, "(",d,")", "RolloutEnter  s=",s, " -> sp =",sp) 
        a  = spw.pars.rolloutPolicy(s,spw.pars.rng)
        #This runs a roll-out simulation, note the lack of discount factor...
        
        R = spw.pars.getReward(s,a, spw.pars)::Reward
        if !spw.pars.terminate
            spw.pars.getNextState!(sp,s,a,spw.pars.rng)
            #println(tabs, "(",d,")", "RolloutNext  s=",s, " -> sp =",sp)
            R += spw.pars.γ*rollout!(spw,sp,s,int16(d-1))::Reward
        end
        
        #println(tabs, "(",d,")", "RolloutExit  s=",s, " -> sp =",sp)
    end 
    return R
end

end # module