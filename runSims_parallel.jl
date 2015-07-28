using pattern
using sim3d


const __POLICY__ = :KRON # :KRON


#Begin by silent policy
loadPolicy = beta -> (S::SType -> pattern.g_noaction)
if __POLICY__ == :MCTS
    using CTMDP_mcts
    loadPolicy = loadMCTSPolicy
elseif __POLICY__ == :KRON
    using CTMDP_kron
    loadPolicy = beta -> loadCTMDPpolicy(1.0, beta)
else
    println("Using silent policy!")
end


function runBatchSimsParallel(seedVal::Int64)
    
    betaVals = [0.0f0, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
    Nbatch = 10 #10?
    tBatchHours = 20.

    if __POLICY__ == :MCTS
       betaVals = [0.0f0, 0.003f0]
       Nbatch = 1 #10?
       tBatchHours = 20.
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

    pars = ["α", pattern.α, "nPhases", pattern.nPhases]
    
    if __POLICY__ == :MCTS 
        append!(pars, ["policy","MCTS", 
                "d", mcts.pars.d, "n", mcts.pars.n,
                "ec", mcts.pars.ec,
                "γ", mcts.pars.γ, "resetDict", mcts.pars.resetDict])
    elseif __POLICY__ == :KRON
        append!(pars, ["policy","KRON"]) 
    end
   return pars                 
end

function fileprefix()
    prefix = string(__POLICY__) * "_n" * string(pattern.nPhases)
    return prefix
end


