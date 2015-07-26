# parallel test code
using HDF5, JLD
using Dates

tstart = now()

parallel = :both

ncpu_local = 0 #cambridge too slow?
machines = [("zouhair@cheonan.stanford.edu", 5, "/usr/bin"),
            ("zouhair@tula.stanford.edu", 5, "/usr/bin")]


#TODO: Figure out how to make this more robust?
sshflags = `-i /home/zouhair/.ssh/id_rsa_cambridge`

if parallel == :local || parallel == :both
    println("Adding ", ncpu_local, " local CPUs")
    addprocs(int64(ncpu_local))
end

if parallel == :remote || parallel == :both
    for (machine, count, dir) in machines
        println("Adding ", count, " " ,machine, " CPUs")
        cluster_list = ASCIIString[]

        for i = 1:count
            push!(cluster_list, machine)
        end

        addprocs(cluster_list, dir = dir, sshflags=sshflags)
    end
end


@everywhere __PARALLEL__ = true

#################################
println("Loading code everywhere")
#################################
require("runSims_parallel.jl")

#################################
println("running using #", pattern.nPhases, " phases")
#################################
tic()
allResults = runAllSims()
concResults = sim3d.concatenate(allResults)

#################################
filename = "mcts_simResults_n"*string(pattern.nPhases)*"_"*string(tstart)*".jld"
println("Saving results to "*filename)
#################################
#Just save all of the data, we'll deal with concatennating later...
JLD.save(filename,  "betaVals", concResults.betaVals,
                    "alertCounts", concResults.alertCounts,
                    "flightTimes", concResults.flightTimes, 
                    "tTotals", concResults.tTotals, 
                    "nNMACcounts", concResults.nNMACcounts
                    "d", mcts.pars.d, "n", mcts.pars.n, "ec", mcts.pars.ec,
                    "γ", mcts.pars.γ,
                    "resetDict", mcts.pars.resetDict)

println("Done !")
toc()    
