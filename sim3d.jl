module sim3d

using airplaneType
using SASS_sensor
using pattern
using posType

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
      if act != pattern.g_noaction && acList[act[1]].readyForATC
        acList[act[1]].atcCommand = act[2]
        alertCount += 1
      end
    end
    
#     if (t > 4640. && readyForCommand)
#         S = [(pattern.appendPhase(ac.navDest[1],ac.navPhase) , ac.readyForATC) for ac in acList]
#         println(t, "", S, "=> \t\t", act)
#         airplaneType.simPars.debugFlag = true
#     end

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



function rstart_in()
    #TODO: find a better way to do this
#     states = pattern.allstates[rand(3:length(pattern.allstates), 4)
    m = 3; M = pattern.g_nNodes;
    idx = int64(rand(pattern.rng, 4) * (M-m) + m)
    return [airplane(46 + (rand(pattern.rng)-0.5)*3,pattern.allstates[i]) for i in idx];
    

end
#################################################
function randomStart()
#TODO: Make this start not at just the nodes!
#################################################
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
    runTime::Float64
    #Outputs
    tTotals    ::Array{Float32, 3}
    flightTimes::Array{Float32, 3}
    alertCounts::Array{Uint32, 3}
    nNMACcounts::Array{Uint32, 3}
    
    collisionPos::Array{Vector{pos}, 3}
    seedValues::Array{Vector{Uint32}, 3}
    
    function simResults(betaVals::Vector{Float32}, atcTypes::Vector{Symbol}, Nbatch)
        NbetaVals = length(betaVals)
        NatcTypes = length(atcTypes)
        new(
        copy(betaVals), copy(atcTypes),
        0., 
        zeros(Float32, NbetaVals,NatcTypes,Nbatch),
        zeros(Float32, NbetaVals,NatcTypes,Nbatch),
        zeros(Uint32 , NbetaVals,NatcTypes,Nbatch),
        zeros(Uint32 , NbetaVals,NatcTypes,Nbatch),
        Array(Vector{pos},NbetaVals, NatcTypes, Nbatch),
        Array(Vector{Uint32}, NbetaVals, NatcTypes, Nbatch)
        )
    end
    
    function simResults(NbetaVals::Number, NatcTypes::Number, Nbatch)
        betaVals = zeros(Float32, NbetaVals)
        atcTypes = fill(:a, NatcTypes)
        return simResults(betaVals, atcTypes, Nbatch)
    end

    
    
end



function concatenate(results)
    res1 = results[1]

    Nprocs = length(results)
    Nbatch = size(res1.flightTimes,3)
    
    concResults = simResults(length(res1.betaVals), length(res1.atcTypes), Nbatch*Nprocs)
    
    concResults.betaVals = res1.betaVals
    concResults.atcTypes = res1.atcTypes
    concResults.runTime = sum([res.runTime for res in results])
    
    cnt = 1:Nbatch
    for i in 1:Nprocs                    
        for field in [:flightTimes, :alertCounts, :nNMACcounts, :tTotals]        
            concResults.(field)[:,:, cnt] = results[i].(field)[:,:,:]            
        end
        cnt += Nbatch
    end
    return concResults
end

#################################################
function runBatchSims(betaVals::Vector{Float32}, 
                      tBatchTime_hours::Number, Nbatch::Int64,
                      seedStart, loadPolicy::Function; 
                      atcTypes::Vector{Symbol} = [:Smart], 
                      Verbosity::Symbol = :None,
                      saveCollPos = false,
                      saveSeeds = false,
                      seedStride = 1000)
#################################################
    # betaVals = [0.0f0, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
    # betaVals = [0.0f0, 0.001f0, 0.005f0] #, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
    # betaVals = [0.0f0]
    
    #atcTypes = [:Smart, :None] , :None is no ATCtype!
    
    seedVal0 = uint32(seedStart * seedStride);
    seedWarning =  uint32(seedVal0 + seedStride);
    
    if Verbosity != :None
        println("Running with seed: ", seedVal0)
    end

    #Total run time per batch in seconds
    tBatchTime = 3600.*tBatchTime_hours
    
    
    #Allocate the result vectors
    results = simResults(betaVals, atcTypes, Nbatch)
    NbetaVals = length(betaVals)
    NatcTypes = length(atcTypes)
    #These should be references just to make hte code a bit more readable
    tTotals      = results.tTotals
    flightTimes  = results.flightTimes
    alertCounts  = results.alertCounts
    nNMACcounts  = results.nNMACcounts
    collisionPos = results.collisionPos
    seedValues   = results.seedValues
    
    
    startTime = time();
    for (betaIdx, beta) in enumerate(betaVals)
        
        #TODO: Consider moving this to inside of the for(atc) loop?
        seedVal = seedVal0; 
        
        
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
                collisionPos[betaIdx,atcIdx, :] = collisionPos[1,atcIdx, :]
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
                                    " -> Perf = $Perf (hr_to_nmac $(round(ft,2)) | $nNmac) \t", 
                                    " Batch = $i / $Nbatch", 
                                    " betaIdx= $betaIdx / $NbetaVals",
                                    " atcIdx= $atcIdx / $NatcTypes")                            
                        end 
                    end   
                    
                    
                    #Reinit seeds for repeatibility
                    if(seedVal == seedWarning)
                        println("seed value has wrapped around the seed stride")
                    end
                    
                    srand(pattern.rng, seedVal)                    
                    if saveSeeds
                        if !isdefined(seedValues,betaIdx,atcIdx, i)
                            seedValues[betaIdx,atcIdx, i] = [seedVal]
                        else
                            push!(seedValues[betaIdx,atcIdx, i], seedVal)
                        end
                    end
                    seedVal += 1

    
                    #First, initialize the position
                    aircraftList = randomStart()
                    #Simulate the policy
                    simTime = tBatchTime-tTotals[betaIdx, atcIdx, i]
                    (idmin, tmax, alertCount, flightTime) = simulate!(aircraftList,simTime,
                                                                      atcType,policy,
                                                                      stopEarly=true, savepath=false, 
                                                                      savemeasurements=false);
                    
                    #we finished early, it must have been an NMAC
                    if (tmax < simTime )
                        nNMACcounts[betaIdx, atcIdx, i] += 1
                    end
                    
                    tTotals[betaIdx, atcIdx, i]     += tmax
                    alertCounts[betaIdx, atcIdx, i] += alertCount
                    flightTimes[betaIdx, atcIdx, i] += flightTime
                    
                    #Keep track of collision loacations
                    if saveCollPos && !(0 in idmin)
                        collPos = copy(aircraftList[idmin[1]].posNED)
                        if !isdefined(collisionPos,betaIdx,atcIdx, i)
                            collisionPos[betaIdx,atcIdx, i] = [collPos]
                        else
                            push!(collisionPos[betaIdx,atcIdx, i], collPos)
                        end
                    end
                    
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


