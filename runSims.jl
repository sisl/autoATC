include("pattern.jl");
include("CTMDP.jl")
include("3dsim.jl"); 

println("Running with nPhases=", nPhases)

α = 1.0;

Aopt = loadPolicy(1.,0.0f0)
ctmdpPolicy = S -> policy_S2a(S, Aopt)


srand(rng,uint32(1));



data_2phase = load("policies/steadystatedistribution_n_2.jld")
const α_π = data_2phase["α_π"]
const g_XDIMS_2phase = data_2phase["g_XDIMS"]
const g_allstates_2phase = data_2phase["g_allstates"]
function LIDX2X_2phase(lidx::Int64)
  #TODO: make this faster ...
  X_rev = xType[ind2sub(g_XDIMS_2phase, lidx) ...]
  return reverse(X_rev)
end
function LIDX2S_2phase(lidx::Int64)
  X = LIDX2X_2phase(lidx)
  return sType[(g_allstates_2phase::Vector{Symbol})[x] for x in X]
end


#rstart_in = ()-> [airplane(46 + (rand()-0.5)*3,s) for s in allstates[rand(3:length(allstates), 4)]];
function rstart_in()
     #we use the 2phase version since the weighted choices come from that ...
     S = LIDX2S_2phase(weightedChoice(α_π))
     return [airplane(46 + (rand(rng)-0.5)*3, s,  (phaseNum(s)-1) / nPhases) for s in S]
     
     #Uniform, ignores the taxi/runway states
     #TODO: fix this
     #return [airplane(46 + (rand()-0.5)*3,s) for s in allstates[rand(3:length(allstates), 4)]];
end

function randomStart()
    
    acList = rstart_in()

    idmin = [0,0]
    notok = (ACLIST)-> getDmin!(idmin, ACLIST)[1] <=  0
        
    rstrt = 0
    while(notok(acList))
        rstrt += 1
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
                                simTime = tBatchTime-tTotals[betaIdx, atcIdx, i]
                (idmin, tmax, alertCount, flightTime) = simulate!(aircraftList,simTime,
                                                                  atcType,ctmdpPolicy,
                                                                  stopEarly=true, savepath=false);
                
                
                if (tmax < simTime ) #we finished early, must have been an NMAC
                    nNMACcounts[betaIdx, atcIdx, i] += 1
                end
                
                tTotals[betaIdx, atcIdx, i] += tmax
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


@time save("./results/steadystate_3DSimResults_n"*string(nPhases)*".jld",  "tStat", tStat, "alertCounts", alertCounts,  "flightTimes", flightTimes, "tTotals", tTotals, "nNMACcounts", nNMACcounts);


