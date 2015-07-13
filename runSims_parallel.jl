using pattern
using CTMDP_mcts
using sim3d


using HDF5, JLD 


function runBatchSimsParallel(seedVal::Int64)
    
    betaVals = [0.0f0, 0.001f0, 0.005f0, 0.01f0]
    Nbatch = 4 #10?
    tBatchHours = 0.1 # 20h?
    return runBatchSims(betaVals, tBatchHours, Nbatch,
                        seedVal, loadMCTSPolicy)

end


function runAllSims()

    #Each processor gets a different seedvalue
    seedVals = [1:(nprocs()-1)]
    
    #Run everything in parallel...
    allResults = pmap(runBatchSimsParallel, seedVals)
    
    #Just save all of the data, we'll deal with concatennating later...
    JLD.save("mcts_simResults.jld", "allResults", allResults)
    
end