using pattern
using sim3d


#const __POLICY__ = :KRON 
#const __POLICY__ = :MCTS
const __POLICY__  = :GSMDP
#const __POLICY__ = :SILENT

#Begin by silent policy
loadPolicy = beta -> ((S::SType, E::Vector{Float32}) -> pattern.g_noaction)

if __POLICY__ == :MCTS
    using CTMDP_mcts
    loadPolicy = loadMCTSPolicy
elseif __POLICY__ == :GSMDP
    using GSMDP_mcts
    loadPolicy = loadMCTSPolicy_gsmdp
elseif __POLICY__ == :KRON
    using CTMDP_kron
    loadPolicy = beta -> loadCTMDPpolicy(1.0, beta)
elseif __POLICY__ == :SILENT
    println("Using silent policy!")
else
    error("Unknown policy ", __POLICY__)
end


function runBatchSimsParallel(seedVal::Int64)
    
    betaVals = [0.0f0, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
    Nbatch = 10 #10?
    tBatchHours = 24.

    if __POLICY__ == :MCTS
       betaVals = [0.0f0, 0.005f0, 0.01f0] # , 0.003f0]
       Nbatch = 1 #10?
       tBatchHours = 24.
    elseif __POLICY__ == :GSMDP
       betaVals = [0.0f0] # , 0.003f0]
       Nbatch = 1 #10?
       tBatchHours = 24.
    elseif __POLICY__ == :SILENT
       betaVals = [0.0f0] # , 0.003f0]
    end

    return runBatchSims(betaVals, tBatchHours, Nbatch,
                        seedVal, loadPolicy; Verbosity=:High)

end


function runAllSims()

    #Each processor gets a different seedvalue
    seedVals = [1:(nprocs()-1)]
    
    #Run everything in parallel...
    allResults = pmap(runBatchSimsParallel, seedVals)
    
    concResults = sim3d.concatenate(allResults)

    return concResults

end

function runPars() 

    pars = ["α", pattern.α, "nPhases", pattern.nPhases, "policy", string(__POLICY__)]
    
    if __POLICY__ == :MCTS 
        append!(pars, [
                "d", mcts.pars.d, "n", mcts.pars.n,
                "ec", mcts.pars.ec,
                "ζ", mcts.pars.ζ, "resetDict", mcts.pars.resetDict])
    end
   return pars                 
end

function fileprefix()
    prefix = string(__POLICY__) * "_n" * string(pattern.nPhases)
    return prefix
end


