# import Base.print
# import Base.show
# import Base.string

using airplaneType

const simdt = 0.25

#################################################
function runAutoATC(acList::Vector{airplane}, policyFun)
#################################################
  act = g_noaction
  S = Symbol[appendPhase(ac.navDest[1],ac.navPhase) for ac in acList]
  act = policyFun(S)
  return act
end


#################################################
function simulate!(acList::Vector{airplane}, Tend, policyTiming::Symbol, policyFun; stopEarly = false, savepath = true)
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
  tidx = 0
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
      if act != g_noaction && acList[act[1]].atcCommand == :∅
        acList[act[1]].atcCommand = act[2]
        alertCount += 1
      end
    end

    #Fly pattern logic for all aircraft
    for idx in 1:length(acList)
      flyPattern!(acList[idx], simdt, savepath)
    end

    #Compute the distance to all other boogies
    dmin = getDmin!(idmin, acList)
    #We had an NMAC event if dmin <= 0, so break out
    if(stopEarly && dmin <= 0)
      break;
    end

    #Accumulate the amount of time spent in non-taxi states
    for idx in 1:length(acList)
      if acList[idx].navDest[1] != :T
        flightTime += simdt
      end
    end

  end

  tmax = Inf
  if(tidx != -1)
    tmax = trange[tidx]
  end
  return (idmin, tmax, alertCount, flightTime/length(acList))
end




