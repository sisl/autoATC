include("pattern.jl");
include("CTMDP.jl")
include("3dsim.jl"); 

println("Running with nPhases=", nPhases)

α = 1.0;

Aopt = loadPolicy(1.,0.0f0)
ctmdpPolicy = S -> policy_S2a(S, Aopt)


srand(rng,uint32(1));

function randomStart()
    rstart_in = ()-> [airplane(46 + (rand()-0.5)*3,s) for s in allstates[rand(3:length(allstates), 4)]];
    acList = rstart_in()

    idmin = [0,0]
    notok = (ACLIST)-> getDmin!(idmin, ACLIST)[1] <=  0
        
    while(notok(acList))
        acList = rstart_in()
    end
    return acList
end


tic()
tBatchTime = 3600*20 #Total run time per batch
N = 100 #Total Number of batches (total simulated time is tBatch * N )

betaVals = [0.0f0, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
betaValStr = [string(β) for β in betaVals]
NbetaVals = length(betaVals)
atcTypes = [:Smart, :None] #Get rid of periodic?
NatcTypes = length(atcTypes)

tStat = zeros(Float32, NbetaVals,NatcTypes,N)
tTotals = zeros(tStat)
flightTimes = zeros(tStat)

alertCounts = zeros(Uint32, NbetaVals,NatcTypes,N)
nNMACcounts = zeros(alertCounts)
toc()
#Last 3 is for 3D position @ collision
#collisionPos = (Symbol => Array{Float64,3})[a => fill(NaN, NbetaVals,N,3) for a in atcTypes]; 


#Reinit sim seed for repeatibility
srand(rng,uint32(1))
startTime = time();
@time for (betaIdx, beta) in enumerate(betaVals)
    #Use the appropriate policy
    s0 = time()
    Aopt = loadPolicy(1., beta)
    @printf("Loading policy took %.0f(s)\n",time()-s0)
    for (atcIdx, atcType) in enumerate(atcTypes)  
        #We only need to run the first beta for the None case
        #Copy it over for simplicity
        if (atcType == :None && betaIdx > 1)
            tStat[betaIdx, atcIdx, :] = tStat[1, atcIdx, :]
            flightTimes[betaIdx, atcIdx, :] = flightTimes[1, atcIdx, :]
            alertCounts[betaIdx, atcIdx, :] = alertCounts[1, atcIdx, :]
            nNMACcounts[betaIdx, atcIdx, :] = nNMACcounts[1, atcIdx, :] 
            tTotals[betaIdx, atcIdx, :] = tTotals[1, atcIdx, :] 

            #collisionPos[atcType][betaIdx,:,:] = collisionPos[atcType][1,:,:]
            continue
        end
        #Simulate N trajectories
        nextPrintTime = time() + 60 
        for i in 1:N
            now = time()
            if now >= nextPrintTime
                nextPrintTime = now + 60
                println(now-startTime, "(s) at betaIdx=",betaIdx, " atcIdx=",atcIdx, " i = ", i)
            end
            
            while tTotals[betaIdx, atcIdx, i] < tBatchTime
                #First, initialize the position
                aircraftList = randomStart()
                #Simulate the policy
                (idmin, tmax, alertCount, flightTime) = simulate!(aircraftList,tBatchTime,atcType,ctmdpPolicy,
                                                                  stopEarly=true, savepath=false);
                tTotals[betaIdx, atcIdx, i] += tmax
                if (tmax < tBatchTime )
                    nNMACcounts[betaIdx, atcIdx, i] += 1
                end
                alertCounts[betaIdx, atcIdx, i] += alertCount
                flightTimes[betaIdx, atcIdx, i] += flightTime
                tStat[betaIdx, atcIdx, i] = tmax
                
                #Whine if for some reason we got a collision right
                #after we started... randomStart() should not cause this
                if(tmax == 0)
                    error("Collision right after start. This shouldn't happen!")
                end
            end

            #Keep track of statistics
#             if !(0 in idmin)
#                 collisionPos[atcType][betaIdx,i,:] = 
#                     [aircraftList[idmin[1]].posNED.n, aircraftList[idmin[1]].posNED.e, aircraftList[idmin[1]].posNED.d]
#             end
        end
  end
end


@time save("./results/new_3DSimResults_n"*string(nPhases)*".jld",  "tStat", tStat, "alertCounts", alertCounts,  "flightTimes", flightTimes, "tTotals", tTotals, "nNMACcounts", nNMACcounts);


