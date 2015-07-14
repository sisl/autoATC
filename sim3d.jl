module sim3d

using airplaneType
using SASS_sensor
using pattern

export simulate!, randomStart
export simResults, runBatchSims

#################################################
#Simulation parameters
#################################################
#This seems way too fine for the resolution we need?
#Would be nice to make this closer to 1second?
const simdt = 0.25

#################################################
function runAutoATC(acList::Vector{airplane}, policyFun)
#################################################
  S = Symbol[pattern.appendPhase(ac.navDest[1],ac.navPhase) for ac in acList]
  return policyFun(S)
end


#################################################
function simulate!(acList::Vector{airplane}, Tend, policyTiming::Symbol, policyFun; 
                   stopEarly = false, savepath = true, savemeasurements = false)
#################################################
#Running simulation,

  assert(policyTiming in [:Smart, :Periodic, :None])

  #Total time range
  trange =  0:simdt:Tend


  stopsim = false;
  idmin = Int64[0,0]
  dmin = Inf
  alertCount = 0
  flightTime = 0.
  
  measurements = []
  
  if(savemeasurements)
    measurements = [ne_psi() for i in 1:length(acList), j in 1:(length(trange)+1)]
    for acIdx in 1:length(acList)
        SASS_sense!(measurements[acIdx,1], acList[acIdx])
    end
  end
  
  tmax = Tend
  for (tidx, t) in enumerate(trange)
    #Find out if any of the aircraft in the pattern
    #is ready for a command. This could be done more
    #concisely but list comprehensions seem to slow
    #things down!
    readyForCommand = false
    if(policyTiming == :Smart)
      #Find out if any aircraft is ready for a command
      for idx in 1:length(acList)
        readyForCommand = readyForCommand || acList[idx].readyForATC
      end
    elseif policyTiming == :Periodic
      #Do it based on clock, once every 10 seconds
      readyForCommand = (t % 10 == 0)
    #else, keep readyForCommand = false, i.e. we won't issue any commands
    end

    #If any aircraft is about to transition,
    #see if there's an ATC command that should be issued!
    if(readyForCommand)
      act = runAutoATC(acList, policyFun)
      #If we have an action to issue, pass it along
      if act != pattern.g_noaction && acList[act[1]].atcCommand == :∅
        acList[act[1]].atcCommand = act[2]
        alertCount += 1
      end
       #println(t, "Ready for command:", [(pattern.appendPhase(ac.navDest[1],ac.navPhase) , ac.readyForATC) for ac in acList], act)
    end

    #Fly pattern logic for all aircraft
    for idx in 1:length(acList)
      flyPattern!(acList[idx], simdt, savepath)
      if savemeasurements
        SASS_sense!(measurements[idx, tidx+1], acList[idx])
      end
    end

    #Compute the distance to all other boogies
    dmin = getDmin!(idmin, acList)
    #We had an NMAC event if dmin <= 0, so break out
    if(stopEarly && dmin <= 0)
      tmax = t;
      break;
    end

    #Accumulate the amount of time spent in non-taxi states
    for idx in 1:length(acList)
      if acList[idx].navDest[1] != :T
        flightTime += simdt
      end
    end

  end

  return (idmin, tmax, alertCount, flightTime/length(acList), measurements)
end

#################################################
function randomStart()
#################################################
    rstart_in = ()-> [airplane(46 + (rand()-0.5)*3,s) for s in pattern.allstates[rand(3:length(pattern.allstates), 4)]];
    acList = rstart_in()

    idmin = [0,0]
    notok = (ACLIST)-> getDmin!(idmin, ACLIST)[1] <=  0
        
    while(notok(acList))
        acList = rstart_in()
    end
    return acList
end



#################################################
#Results from running batch sims
type simResults
    #Inputs
    betaVals::Vector{Float32}
    atcTypes::Vector{Symbol}
    seedVal::Int64
    runTime::Float64
    #Outputs
    tTotals    ::Array{Float32, 3}
    flightTimes::Array{Float32, 3}
    alertCounts::Array{Uint32, 3}
    nNMACcounts::Array{Uint32, 3}
    
    
    function simResults(betaVals, atcTypes, Nbatch, seedVal)
        NbetaVals = length(betaVals)
        NatcTypes = length(atcTypes)
        new(
        copy(betaVals), copy(atcTypes),
        seedVal, 0., 
        zeros(Float32, NbetaVals,NatcTypes,Nbatch),
        zeros(Float32, NbetaVals,NatcTypes,Nbatch),
        zeros(Uint32 , NbetaVals,NatcTypes,Nbatch),
        zeros(Uint32 , NbetaVals,NatcTypes,Nbatch)
        )
    end
    
    function simResults(NbetaVals, NatcTypes, Nbatch)
        new(
        zeros(Float32, NbetaVals),
        Array(Symbol , NatcTypes),
        0, 0., 
        zeros(Float32, NbetaVals,NatcTypes,Nbatch),
        zeros(Float32, NbetaVals,NatcTypes,Nbatch),
        zeros(Uint32 , NbetaVals,NatcTypes,Nbatch),
        zeros(Uint32 , NbetaVals,NatcTypes,Nbatch)
        )
    end

    
    
end


function concatenate(results)
    res1 = results[1]

    Nprocs = length(results)
    concResults = simResults(length(res1.betaVals), length(res1.atcTypes), size(res1.flightTimes, 3)*Nprocs)
    
    concResults.betaVals = res1.betaVals
    concResults.atcTypes = res1.atcTypes
    concResults.runTime = sum([res.runTime for res in allResults])
    
    for field in [:flightTimes, :alertCounts, :nNMACcounts, :tTotals]
        cnt = 1:2
        for i in 1:Nprocs
            assert(res1.betaVals == results[i].betaVals)
            assert(res1.atcTypes == results[i].atcTypes)
            concResults.(field)[:,:, cnt] = results[i].(field)[:,:,:]
            cnt += 2
        end
    end
    return concResults
end

#################################################
function runBatchSims(betaVals::Vector{Float32}, 
                      tBatchTime_hours::Number, Nbatch::Int64,
                      seedVal::Int64, loadPolicy::Function; 
                      atcTypes::Vector{Symbol} = [:Smart], Verbosity::Symbol = :None)
#################################################
    # betaVals = [0.0f0, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
    # betaVals = [0.0f0, 0.001f0, 0.005f0] #, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
    # betaVals = [0.0f0]
    
    #atcTypes = [:Smart, :None] , :None is no ATCtype!
    
    if Verbosity != :None
        println("Running with seed: ", seedVal)
    end

    #Reinit seeds for repeatibility
    #Note that although we are using the same seedValue for these RNGs,
    #they are used in different contexts so it doesn't really matter
    srand(pattern.rng, uint32(seedVal))
    #srand(CTMDP_mcts.mctsRng, uint32(seedVal))
    
    #Total run time per batch in seconds
    tBatchTime = 3600.*tBatchTime_hours
    
    
    #Allocate the result vectors
    results = simResults(betaVals, atcTypes, Nbatch, seedVal)
    NbetaVals = length(betaVals)
    NatcTypes = length(atcTypes)
    #These should be references just to make hte code a bit more readable
    tTotals     = results.tTotals
    flightTimes = results.flightTimes
    alertCounts = results.alertCounts
    nNMACcounts = results.nNMACcounts
    
    
    startTime = time();
    for (betaIdx, beta) in enumerate(betaVals)
        #Use the appropriate policy
        s0 = time()
        policy = loadPolicy(beta)::Function
        for (atcIdx, atcType) in enumerate(atcTypes)  
            #We only need to run the first beta for the None case
            #and copy it over for simplicity
            if (atcType == :None && betaIdx > 1)
                flightTimes[betaIdx, atcIdx, :] = flightTimes[1, atcIdx, :]
                alertCounts[betaIdx, atcIdx, :] = alertCounts[1, atcIdx, :]
                nNMACcounts[betaIdx, atcIdx, :] = nNMACcounts[1, atcIdx, :] 
                tTotals[betaIdx, atcIdx, :] = tTotals[1, atcIdx, :] 
    
                #collisionPos[atcType][betaIdx,:,:] = collisionPos[atcType][1,:,:]
                continue
            end
            nextPrintTime = time() 
            
            #Simulate Nbatch_es
            for i in 1:Nbatch
                
                #Each batch should run all the way to tBatchTime!
                while tTotals[betaIdx, atcIdx, i] < tBatchTime
                    if Verbosity != :None
                        now = time()
                        if now >= nextPrintTime
                            nextPrintTime = now + 60
                            
                            ft = sum(flightTimes[betaIdx, atcIdx,1:i])/3600
                            nNmac = sum(nNMACcounts[betaIdx, atcIdx,1:i])
                            Perf = round(ft/nNmac,2)
                            
                            sim_so_far = round(tTotals[betaIdx, atcIdx, i]/3600, 2 )
                            total_to_end = round(Nbatch*tBatchTime/3600, 2)
                            t = round((now-startTime)/60,2)
                            println(" t = $t (mins) ", 
                                    " Simulated Hours = $sim_so_far  / $total_to_end",
                                    " -> Perf = $Perf (hr_to_nmac $(round(ft,2)) | $) \t", 
                                    " Batch = $i / $Nbatch", 
                                    " betaIdx= $betaIdx / $NbetaVals",
                                    " atcIdx= $atcIdx / $NatcTypes")                            
                        end 
                    end   
                    
                    #First, initialize the position
                    aircraftList = randomStart()
                    #Simulate the policy
                    simTime = tBatchTime-tTotals[betaIdx, atcIdx, i]
                    (idmin, tmax, alertCount, flightTime) = simulate!(aircraftList,simTime,
                                                                      atcType,policy,
                                                                      stopEarly=true, savepath=false, savemeasurements=false);
                    
                    #we finished early, it must have been an NMAC
                    if (tmax < simTime )
                        nNMACcounts[betaIdx, atcIdx, i] += 1
                    end
                    
                    tTotals[betaIdx, atcIdx, i]     += tmax
                    alertCounts[betaIdx, atcIdx, i] += alertCount
                    flightTimes[betaIdx, atcIdx, i] += flightTime
                    
                    #Whine if for some reason we got a collision right
                    #after we started... randomStart() should not cause this
                    if(tmax == 0)
                        error("Collision right after start. This shouldn't happen!")
                    end
                end
    
            end
      end
    end
    
    results.runTime = time() - startTime
    
    if Verbosity != :None
        println("Finished processing (total runtime: $(round(results.runTime/60,2)) minutes)")
    end
    
    return results

end



end


