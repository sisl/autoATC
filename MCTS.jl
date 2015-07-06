module MCTS

# This module implements Monte Carlo tree search (MCTS), an online MDP solver. In MCTS. the poosible evolutions of the system are represented as a tree and an approximation to this tree is built iteratively and used to inform the choice of action. Documentation for MCTS can be found in Decision Making Under Uncertainty, MIT Press 2015, and A survey of Monte Carlo tree search methods in Computational Intelligence and AI in games, 2012.

using MDP
using auxfuncs

export SPWParams, SPW, selectAction

typealias Depth Int16

type SPWParams{T<:Action}
    d::Depth                    # search depth
    ec::Float64                 # exploration constant- governs trade-off between exploration and exploitation in MCTS
    n::Int32                    # number of iterations
    rng::AbstractRNG            # random number generator
    A::Array{T,1}               # set of allowable actions
    getAction::Function         # returns action for rollout policy
    getNextState::Function      # takes state and action as arguments and returns next state from generative model
    getReward::Function         # takes state and action as arguments and returns reward
end

type StateNode
    n::Array{Int32,1}
    q::Array{Reward,1}
    StateNode(nA) = new(zeros(Int32,nA),zeros(Reward,nA))
end

type SPW{T}
    s::Dict{State,StateNode}
    p::SPWParams{T}
    SPW{T<:Action}(p::SPWParams{T}) = new(Dict{State,StateNode}(),p)
end

function selectAction(spw::SPW,s::State)
    # This function calls simulate and chooses the approximate best action from the reward approximations 
    for i = 1:spw.p.n 
        simulate(spw,s,spw.p.d)
    end
    return spw.p.A[indmax(spw.s[s].q)]::Action # Choose action with highest apporoximate value
end

function simulate(spw::SPW,s::State,d::Depth)
    # This function returns the reward for one iteration of MCTS
    if d == 0
        return 0.0::Reward
    end
    if !haskey(spw.s,s)
        spw.s[s] = StateNode(length(spw.p.A))
        return rollout(spw,s,d)::Reward
    else # choose an action using UCT criterion
        cS = spw.s[s]
        i = indmax(cS.q + spw.p.ec.*real(sqrt(complex(log(sum(cS.n))./cS.n))))
        a = spw.p.A[i] # choose action with highest UCT score
        sp = spw.p.getNextState(s,a,spw.p.rng)
        q = spw.p.getReward(s,a) + simulate(spw,sp,int16(d-1))
        cS.n[i] = cS.n[i] + one(Int32)
        cS.q[i] += (q-cS.q[i])/cS.n[i]
        return q::Reward
    end
end

function rollout(spw::SPW,s::State,d::Depth)
    # Runs a rollout simulation using the default policy
    if d == 0
        return 0.0::Reward
    else 
        a = spw.p.getAction(s,spw.p.rng)
        sp = spw.p.getNextState(s,a,spw.p.rng)
        return (spw.p.getReward(s,a) + rollout(spw,sp,int16(d-1)))::Reward
    end 
end

end # module