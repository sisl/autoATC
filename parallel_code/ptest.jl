# parallel test code

parallel = :both
ncpu_local = CPU_CORES / 2
machines = [("zouhair@cheonan.stanford.edu", 4, "/usr/bin"), ("zouhair@tula.stanford.edu", 20, "/usr/bin")]

if parallel == :local_ || parallel == :both
    addprocs(int64(ncpu_local))
end

if parallel == :remote || parallel == :both
    for (machine, count, dir) in machines
        cluster_list = ASCIIString[]

        for i = 1:count
            push!(cluster_list, machine)
        end

        addprocs(cluster_list, dir = dir)
    end
end


@everywhere __PARALLEL__ = true
require("pworker.jl")


runTest(32, bParallel = true)


