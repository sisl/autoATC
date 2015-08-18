module MCTS_GSMDP
#This module implements Monte Carlo tree search (MCTS) applied to GSMDP
#It is an adaptation of sparse UCT search to continuous time


export SPWParams, SPW, selectAction!

typealias Depth Int16
typealias Reward Float32


#TODO: This forces a certain type for state and action ...
typealias State       Vector{Symbol}
typealias EventLength Vector{Float32}
typealias Action (Int8, Symbol)
typealias StateKey (Vector{Symbol}, Vector{Float32}, Vector{Int8})
typealias StateEvent StateKey


type SPWParams{T<:Action}
    terminate::Bool             #Flag to be populated by getReward to indicate we have reached terminal state
    resetDict::Bool             #Whether to reset dictionary
    buildTree::Bool             #Whether to build tree for debugging purposes

    d::Depth                    # search depth
    ec::Float32                 # exploration constant- governs trade-off between exploration and exploitation in MCTS
    n::Int32                    # number of iterations
    β::Float32                  # Alert/Collision ration (should be inside of problem defintion...)
    ζ::Float32                  # Discount rate
    
    w::Int64                    # Width for sparse UCT sampling! 
    
    Afun!::Function             # set of allowable actions 
    rolloutPolicy::Function     # returns action for rollout policy
    getNextState!::Function     # takes state and action as arguments and populates next state from generative model
    getReward::Function         # takes state and action as arguments and returns reward
    hashState::Function
    rng::AbstractRNG            # random number generator

end


#statistics for a given node
type StateStat
    n::Vector{Float32} #Number of times the state/action pairs were visitied
    q::Vector{Reward}  #average reward for that state
    
    t_sojurn::Array{Float32,  2}
    children::Array{StateKey, 2}
    childrenCnt::Array{Float32, 2}
    
    StateStat(nA,w) = new(
                            zeros(Float32,nA),
                            zeros(Reward,nA),
                            zeros(Float32,  nA, w),
                            Array(StateKey, nA, w),
                            zeros(Float32,  nA, w))
end 



#Using our own hash to avoid weirdness with dictionaries...
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
function selectAction!(spw::SPW, acts::Vector{Action}, s0t0a0::StateEvent)
       
    nActs = spw.pars.Afun!(acts, s0t0a0)
    
    #If there's only one action to be taken, no point in wasting time...
    if(nActs == 1)
        return acts[1]
    end

    #Reset the dictionary if told to
    if(spw.pars.resetDict)
        #the performance of empty!(spw.stats)
        #is about ~50% better than letting the dictionary grow!
        empty!(spw.stats)
        #spw.stats = Dict{StateKey,StateStat}()
#         for v in values(spw.stats)
#             fill!(v.n, 0.0f0)
#             fill!(v.q, 0.0f0)
#         end
    end    
    
    #This is to avoid the first call to simulate wasting a rollout
    if !haskey(spw.stats, s0t0a0)
        spw.stats[s0t0a0] = StateStat(nActs, spw.pars.w)
    end



    
    for i = 1:spw.pars.n 
        simulate!(spw, acts, s0t0a0 , spw.pars.d)
    end

    #get the list of allowable actions for this state
    #We need to do this since simulate! will work on acts
    #TODO: Find a way to avoid calling this function a 2nd time!
    
    nActs = spw.pars.Afun!(acts, s0t0a0)     
    #Choose action with highest apporoximate q-value 
    return acts[indmax(spw.stats[s0t0a0].q)]::Action 
end



#Sp is passed-in to avoid to have to allocate memory
#Note that this function will mangle s
function simulate!(spw::SPW, acts::Vector{Action}, s_t::StateEvent, d::Depth)
    # This function returns the reward for one iteration of MCTS
    
    #We have reached the bottom, bubble up.
    if d == 0
        return 0.0f0::Reward
    end
    
#      tabs = string([" " for i in 1:(2*d)]...)

    
    #println(tabs, "(",d,")", "Entry s=",s, " -> sp =",sp)
    
    #TODO: Make the A(s) function better!
    #Determine actions available for this state
    nActs = spw.pars.Afun!(acts, s_t)
    
    #nActs = length(acts)::Int64
    
    #If this state has no statistics yet (i.e. first visit)
    #perform a roll-out
    if !haskey(spw.stats,s_t)
        #Add a new node with the state statistics
        spw.stats[s_t] = StateStat(nActs, spw.pars.w)
        #perform one rollout and return
        return rollout!(spw, deepcopy(s_t), deepcopy(s_t), d)::Reward
    else # choose an action using UCT criterion
        # choose action with highest UCT score
        # which trades-off exploration / explotation
        cS = spw.stats[s_t]
        
        N = sum(cS.n)::Float32
        #Check N > 0 to avoid log(0)/0 = -Inf, can't take sqrt!
        iAct = 1 #if N == 0, use first action
        if( N > 0)
            ec_sqrt_logN =  (spw.pars.ec * sqrt(log(N) + eps(0.0f0)))::Float32 # +eps is to avoid log(1)/0 = NaN
            
            #iAct  = indmax(cS.q + ec_sqrt_logN ./ sqrt(cS.n))
            maxQec = -Inf32
            q_ec = 0.0f0
            for j = 1:nActs
                q_ec = (cS.q[j] +  ec_sqrt_logN / sqrt(cS.n[j]))::Reward
                if q_ec > maxQec
                    iAct = j
                    maxQec = q_ec
                end
            end
        end
        a  = acts[iAct]
        
        #Estimate the reward
        q = spw.pars.getReward(s_t, a, spw.pars)::Reward 
                
        #Randomly select the next state based on the action used
        t_sojurn = 0f0;

        if isdefined(cS.children, iAct, spw.pars.w)
            sp_idx = int(floor(rand(spw.pars.rng) * spw.pars.w + 1)) #garanteed <= w
            t_sojurn = cS.t_sojurn[iAct, sp_idx]       
            s_t_p    = deepcopy(cS.children[iAct, sp_idx])
            cS.childrenCnt[iAct, sp_idx] += 1
        else
            s_t_p = deepcopy(s_t)
            t_sojurn = spw.pars.getNextState!(s_t_p, s_t, a,spw.pars.rng)
            for w_i in 1:spw.pars.w
                if !isdefined(cS.children, iAct, w_i) #at least one (spw.pars.w) will not be dined!
                    cS.t_sojurn[iAct, w_i] = t_sojurn
                    cS.children[iAct, w_i] = s_t_p
                    cS.childrenCnt[iAct, w_i] += 1
                    break
                else #check in case we already encountered this!
                    if cS.children[iAct, w_i] == s_t_p
                        cS.childrenCnt[iAct, w_i] += 1
                        break
                    end
                end
            end        
        end
        


#         if(spw.pars.buildTree)
#             cS.STree[i][sp_key] = 1 + get(cS.STree[i], sp_key, 0)
#         end
                
        #println(tabs, "(",d,")", "Next  s=",s, " -> sp =",sp)

        if !spw.pars.terminate
            #Note that a call to simulate! will change s and sp, but we 
            #don't care at this point since we no longer need their values!
            #(s will be used as storage variable for s' in the recursive call)
            q += exp(-spw.pars.ζ*t_sojurn)*simulate!(spw, acts, s_t_p ,int16(d-1))::Reward
            q::Reward
        end
        #println(tabs, "(",d,")", "Exit  s=",s, " -> sp =",sp)

        #Update the statistics
        cS.n[iAct] += 1.0f0
 
        #This could maybe take into account the transition probablities?
        #i.e. given that we know P(s' | s, a), use importance weighing?
        #Doh, apparently people thought about this already. See:
        #http://www.ai.rug.nl/~mwiering/Tom_van_der_Kleij_Thesis.pdf
        #Equ 3.2
        cS.q[iAct] += (q-cS.q[iAct])/cS.n[iAct]
        return q::Reward
    end
end

function rollout!(spw::SPW, s_t::StateEvent, s_t_p::StateEvent, d::Depth)
    # Runs a rollout simulation using the default policy
#      tabs = string([" " for i in 1:(2*d)]...)

    R = 0f0::Reward
    if d > 0
        #println(tabs, "(",d,")", "RolloutEnter  s=",s, " -> sp =",sp) 
        a  = spw.pars.rolloutPolicy(s_t,spw.pars.rng)
        
        R = spw.pars.getReward(s_t, a, spw.pars)::Reward
        if !spw.pars.terminate
            t_sojurn = spw.pars.getNextState!(s_t_p, s_t, a,spw.pars.rng)
            #println(tabs, "(",d,")", "RolloutNext  s=",s, " -> sp =",sp)
            R += exp(-spw.pars.ζ*t_sojurn)*rollout!(spw,s_t_p, s_t, int16(d-1))::Reward
        end
        
        #println(tabs, "(",d,")", "RolloutExit  s=",s, " -> sp =",sp)
    end 
    return R
end

end # module