# parallel test code

parallel = :both
ncpu_local = CPU_CORES / 2
machines = [("zouhair@cheonan.stanford.edu", 4, "/usr/bin"), ("zouhair@tula.stanford.edu", 20, "/usr/bin")]


#TODO: Figure out how to make this more robust?
sshflags = `-i /home/zouhair/.ssh/id_rsa_cambridge`

if parallel == :local_ || parallel == :both
    println("Adding ", ncpu_local, "local CPUs")
    addprocs(int64(ncpu_local))
end

if parallel == :remote || parallel == :both
    for (machine, count, dir) in machines
        println("Adding ", count, machine, " CPUs")
        cluster_list = ASCIIString[]

        for i = 1:count
            push!(cluster_list, machine)
        end

        addprocs(cluster_list, dir = dir, sshflags=sshflags)
    end
end


@everywhere __PARALLEL__ = true

println("Loading code everywere")

require("pworker.jl")

println("running!")
runTest(32, bParallel = true)


