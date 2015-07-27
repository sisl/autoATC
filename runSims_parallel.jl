using pattern
using sim3d



#Begin by silent policy
loadPolicy = beta -> (S::SType -> pattern.g_noaction)
if __POLICY__ == :MCTS
    using CTMDP_mcts
    loadPolicy = loadMCTSPolicy
elseif __POLICY__ == :KRON
    using CTMDP_kronsolver
    loadPolicy = beta -> loadCTMDPpolicy(1.0, beta)
else
    println("Using silent policy!")
end


function runBatchSimsParallel(seedVal::Int64)
    
    betaVals = [0.0f0, 0.001f0]#, 0.005f0, 0.01f0]
    Nbatch = 1 #10?
    tBatchHours = 10
    return runBatchSims(betaVals, tBatchHours, Nbatch,
                        seedVal, ; Verbosity=:High)

end


function runAllSims()

    #Each processor gets a different seedvalue
    seedVals = [1:(nprocs()-1)]
    
    #Run everything in parallel...
    allResults = pmap(runBatchSimsParallel, seedVals)
    
    concResults = sim3d.concatenate(allResults)

    return concResults

end
