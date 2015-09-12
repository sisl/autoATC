# Author: Youngjun Kim, youngjun@stanford.edu
# Date: 12/11/2014

module MCTSVisualizer_

export MCTSVisualizer, initTree, updateTree, saveTree


# using POMDP_

using JSON


type MCTSVisualizer

    b_sim::Bool
    b_hist::Bool
    b_hist_acc::Bool

    sim_file::ASCIIString
    hist_file::ASCIIString
    hist_acc_file::ASCIIString

    T_sim::Dict{ASCIIString, Any}
    T_sim_curr::Dict{ASCIIString, Any}
    T_sim_stack::Vector{Dict{ASCIIString, Any}}

    T_hist::Dict{ASCIIString, Any}
    T_hist_curr::Dict{ASCIIString, Any}
    T_hist_stack::Vector{Dict{ASCIIString, Any}}

    T_hist_acc::Dict{ASCIIString, Any}
    T_hist_acc_curr::Dict{ASCIIString, Any}
    T_hist_acc_stack::Vector{Dict{ASCIIString, Any}}


    function MCTSVisualizer(;sim::Bool = true, sim_file::ASCIIString = "sim.json", hist::Bool = true, hist_file::ASCIIString = "hist.json", hist_acc::Bool = false, hist_acc_file::ASCIIString = "hist_acc.json")

        self = new()

        self.b_sim = sim
        self.b_hist = hist
        self.b_hist_acc = hist_acc

        self.sim_file = sim_file
        self.hist_file = hist_file
        self.hist_acc_file = hist_acc_file

        return self
    end
end


function initTree(vis::MCTSVisualizer)

    if vis.b_sim
        vis.T_sim = Dict{ASCIIString, Any}()
        vis.T_sim_curr = vis.T_sim
        vis.T_sim_stack = Dict{ASCIIString, Any}[]
    end

    if vis.b_hist
        T_hist = Dict{ASCIIString, Any}()
        T_hist["N"] = 0
        T_hist["actions"] = Dict{ASCIIString, Any}()
        vis.T_hist = T_hist
        vis.T_hist_curr = vis.T_hist
        vis.T_hist_stack = Dict{ASCIIString, Any}[]
    end

    if vis.b_hist_acc
        T_hist_acc = Dict{ASCIIString, Any}()
        T_hist_acc["N"] = 0
        T_hist_acc["actions"] = Dict{ASCIIString, Any}()
        vis.T_hist_acc = T_hist_acc
        vis.T_hist_acc_curr = vis.T_hist_acc
        vis.T_hist_acc_stack = Dict{ASCIIString, Any}[]
    end
end


function updateTree(vis::MCTSVisualizer, where::Symbol, args...)

    if where == :start_sim
        if vis.b_sim
            s, = args
            st = string(s)

            T_sim_curr = vis.T_sim_curr

            if !haskey(T_sim_curr, st)
                T_sim_curr[st] = Dict{ASCIIString, Any}()
                T_sim_curr[st]["actions"] = Dict{ASCIIString, Any}()
                T_sim_curr[st]["N"] = 1
            else
                T_sim_curr[st]["N"] += 1
            end
        end

    elseif where == :before_sim
        s, a, o = args

        st = string(s)
        act = string(a.action)
        obs = string(o.observation)

        if vis.b_sim
            T_sim_curr = vis.T_sim_curr

            if !haskey(T_sim_curr[st]["actions"], act)
                T_sim_curr[st]["actions"][act] = Dict{ASCIIString, Any}()
                T_sim_curr[st]["actions"][act]["observations"] = Dict{ASCIIString, Any}()
                T_sim_curr[st]["actions"][act]["N"] = 0
                T_sim_curr[st]["actions"][act]["r"] = 0.
            end

            if !haskey(T_sim_curr[st]["actions"][act]["observations"], obs)
                T_sim_curr[st]["actions"][act]["observations"][obs] = Dict{ASCIIString, Any}()
                T_sim_curr[st]["actions"][act]["observations"][obs]["states"] = Dict{ASCIIString, Any}()
            end

            push!(vis.T_sim_stack, T_sim_curr)
            vis.T_sim_curr = T_sim_curr[st]["actions"][act]["observations"][obs]["states"]
        end

        if vis.b_hist
            T_hist_curr = vis.T_hist_curr

            if !haskey(T_hist_curr["actions"], act)
                T_hist_curr["actions"][act] = Dict{ASCIIString, Any}()
                T_hist_curr["actions"][act]["N"] = 0
                T_hist_curr["actions"][act]["Q"] = 0.
                T_hist_curr["actions"][act]["observations"] = Dict{ASCIIString, Any}()
            end

            if !haskey(T_hist_curr["actions"][act]["observations"], obs)
                T_hist_curr["actions"][act]["observations"][obs] = Dict{ASCIIString, Any}()
                T_hist_curr["actions"][act]["observations"][obs]["N"] = 0
                T_hist_curr["actions"][act]["observations"][obs]["actions"] = Dict{ASCIIString, Any}()
            end

            push!(vis.T_hist_stack, T_hist_curr)
            vis.T_hist_curr = T_hist_curr["actions"][act]["observations"][obs]
        end

        if vis.b_hist_acc
            T_hist_acc_curr = vis.T_hist_acc_curr

            if !haskey(T_hist_acc_curr["actions"], act)
                T_hist_acc_curr["actions"][act] = Dict{ASCIIString, Any}()
                T_hist_acc_curr["actions"][act]["N"] = 0
                T_hist_acc_curr["actions"][act]["Q"] = 0.
                T_hist_acc_curr["actions"][act]["observations"] = Dict{ASCIIString, Any}()
            end

            if !haskey(T_hist_acc_curr["actions"][act]["observations"], obs)
                T_hist_acc_curr["actions"][act]["observations"][obs] = Dict{ASCIIString, Any}()
                T_hist_acc_curr["actions"][act]["observations"][obs]["N"] = 0
                T_hist_acc_curr["actions"][act]["observations"][obs]["actions"] = Dict{ASCIIString, Any}()
            end

            push!(vis.T_hist_acc_stack, T_hist_acc_curr)
            vis.T_hist_acc_curr = T_hist_acc_curr["actions"][act]["observations"][obs]
        end

    elseif where == :after_sim
        s, a, r, q, N, Ns, Q = args

        st = string(s)
        act = string(a.action)

        if vis.b_sim
            vis.T_sim_curr = pop!(vis.T_sim_stack)
            T_sim_curr = vis.T_sim_curr

            T_sim_curr[st]["actions"][act]["N"] += 1
            T_sim_curr[st]["actions"][act]["r"] += (r - T_sim_curr[st]["actions"][act]["r"]) / T_sim_curr[st]["actions"][act]["N"]
        end

        if vis.b_hist
            vis.T_hist_curr = pop!(vis.T_hist_stack)
            T_hist_curr = vis.T_hist_curr

            T_hist_curr["actions"][act]["N"] += 1
            T_hist_curr["N"] += 1
            T_hist_curr["actions"][act]["Q"] += (q - T_hist_curr["actions"][act]["Q"]) / T_hist_curr["actions"][act]["N"]
        end

        # TODO construct T_hist_acc incrementally
        if vis.b_hist_acc
            vis.T_hist_acc_curr = pop!(vis.T_hist_acc_stack)
            T_hist_acc_curr = vis.T_hist_acc_curr

            T_hist_acc_curr["actions"][act]["N"] = N
            T_hist_acc_curr["N"] = Ns
            T_hist_acc_curr["actions"][act]["Q"] = Q
        end
    end
end


function saveTree(vis::MCTSVisualizer, pm::POMDP)

    if vis.b_sim
        saveSimTree(vis, pm)
    end

    if vis.b_hist
        saveHistTree(vis, pm)
    end

    if vis.b_hist_acc
        saveHistTree(vis, pm, acc = true)
    end
end


function saveSimTree(vis::MCTSVisualizer, pm::POMDP)

    function process(Tin, Tout, level; r_prev = 0.)

        if rem(level, 3) == 0
            for (state, node) in Tin
                node_ = Dict{ASCIIString, Any}()
                node_["state"] = state
                node_["N"] = node["N"]
                node_["actions"] = Dict{ASCIIString, Any}[]

                push!(Tout, node_)

                process(node["actions"], node_["actions"], level + 1, r_prev = r_prev)
            end
        elseif rem(level, 3) == 1
            for a in pm.actions
                action = string(a.action)
                if haskey(Tin, action)
                    node = Tin[action]
                else
                    continue
                end

                node_ = Dict{ASCIIString, Any}()
                node_["action"] = action
                node_["N"] = node["N"]
                node_["r"] = node["r"]
                node_["R"] = r_prev + node["r"]
                node_["observations"] = Dict{ASCIIString, Any}[]

                push!(Tout, node_)

                process(node["observations"], node_["observations"], level + 1, r_prev = node_["R"])
            end
        elseif rem(level, 3) == 2
            for o in pm.observations
                observation = string(o.observation)
                if haskey(Tin, observation)
                    node = Tin[observation]
                else
                    continue
                end

                node_ = Dict{ASCIIString, Any}()
                node_["observation"] = observation
                node_["states"] = Dict{ASCIIString, Any}[]

                push!(Tout, node_)

                process(node["states"], node_["states"], level + 1, r_prev = r_prev)
            end
        end
    end

    Tout = Dict{ASCIIString, Any}()
    Tout["name"] = "root"
    Tout["states"] = Dict{ASCIIString, Any}[]

    process(vis.T_sim, Tout["states"], 0)

    f = open(vis.sim_file, "w")
    JSON.print(f, Tout, 2)
    close(f)

    return Tout
end


function saveHistTree(vis::MCTSVisualizer, pm::POMDP; acc::Bool = false)

    function process(Tin, Tout, level)

        if rem(level, 2) == 0
            for a in pm.actions
                action = string(a.action)
                if haskey(Tin["actions"], action)
                    node = Tin["actions"][action]
                else
                    continue
                end

                node_ = Dict{ASCIIString, Any}()
                node_["action"] = action
                node_["N"] = node["N"]
                node_["Q"] = node["Q"]
                node_["observations"] = Dict{ASCIIString, Any}[]

                push!(Tout["actions"], node_)

                process(node, node_, level + 1)
            end
        elseif rem(level, 2) == 1
            for o in pm.observations
                observation = string(o.observation)
                if haskey(Tin["observations"], observation)
                    node = Tin["observations"][observation]

                    if node["N"] == 0
                        continue
                    end
                else
                    continue
                end

                node_ = Dict{ASCIIString, Any}()
                node_["observation"] = observation
                node_["N"] = node["N"]
                node_["actions"] = Dict{ASCIIString, Any}[]

                push!(Tout["observations"], node_)

                process(node, node_, level + 1)
            end
        end
    end

    b_acc = acc

    if !b_acc
        T_hist = vis.T_hist
        hist_file = vis.hist_file
    else
        T_hist = vis.T_hist_acc
        hist_file = vis.hist_acc_file
    end

    Tout = Dict{ASCIIString, Any}()
    Tout["name"] = "root"
    Tout["N"] = T_hist["N"]
    Tout["actions"] = Dict{ASCIIString, Any}[]

    process(T_hist, Tout, 0)

    f = open(hist_file, "w")
    JSON.print(f, Tout, 2)
    close(f)

    return Tout
end


end


