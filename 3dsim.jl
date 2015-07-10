# import Base.print
# import Base.show
# import Base.string

using airplaneType
using SASS_sensor

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
#       println(t, "Ready for command:", [(pattern.appendPhase(ac.navDest[1],ac.navPhase) , ac.readyForATC) for ac in acList], act)
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




