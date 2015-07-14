using pattern
using CTMDP_mcts
using sim3d

function runBatchSimsParallel(seedVal::Int64)
    
    betaVals = [0.0f0]# 0.001f0, 0.005f0, 0.01f0]
    Nbatch = 1 #10?
    tBatchHours = 10 # 5h
    return runBatchSims(betaVals, tBatchHours, Nbatch,
                        seedVal, loadMCTSPolicy)

end


function runAllSims()

    #Each processor gets a different seedvalue
    seedVals = [1:(nprocs()-1)]
    
    #Run everything in parallel...
    allResults = pmap(runBatchSimsParallel, seedVals)
    
    return allResults

end