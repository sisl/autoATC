using HDF5, JLD
using Dates

tstart_string = strftime("%y_%m_%d-%H_%M_%S",time())

parallel = :remote

ncpu_local = 0 #cambridge too slow?

machines = [("zouhair@cheonan.stanford.edu", 8, "/usr/bin"),
            ("zouhair@tula.stanford.edu", 8, "/usr/bin")]


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


#################################
println("Loading code everywhere")
#################################
require("runSims_parallel.jl")

#################################
println("running using #", pattern.nPhases, " phases")
#################################
tic()
concResults = runAllSims()
                    
                    
#################################
filename = "simResults/"*fileprefix()*"_"*tstart_string*".jld"

println("Saving results to "*filename)
#################################
JLD.save(filename,  "betaVals", concResults.betaVals,
                    "alertCounts", concResults.alertCounts,
                    "flightTimes", concResults.flightTimes, 
                    "tTotals", concResults.tTotals, 
                    "nNMACcounts", concResults.nNMACcounts,
                    runPars()...)

println("Done !")
toc()    
