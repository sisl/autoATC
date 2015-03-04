# import Base.print
# import Base.show
# import Base.string

function unwrap(a::Float64)
    if (a > pi)
        return mod(a + pi, 2*pi) - pi
    elseif (a < -pi)
        return -(mod(-a + pi, 2*pi) - pi)  #aka -unwrap(-a)
    else
        return a;
    end
end

function clip(r::Float64, min::Float64, max::Float64)
  if r < min
    return min
  elseif r > max
    return max
  else
    return r
  end
end


#Parameters
const simdt=0.25
const transThresh = 50.
const maxNoise = 150. #FIXME: Make this configurable per airplane? ~500 ft
const maxNoiseAlt = 20. #FIXME: Make this configurable per airplane? ~60 ft

const taxiSpeed = 5. #default taxi speed


#################################################
type pos
#################################################
  n::Float64
  e::Float64
  d::Float64
end

function bearing(p0::pos, p1::pos)
    dN = p1.n - p0.n;
    dE = p1.e - p0.e;
    return unwrap(pi/2-atan2(dN, dE))
end


function distance2(p0::pos, p1::pos)
    dN = p0.n - p1.n;
    dE = p0.e - p1.e;
    return (dN*dN + dE*dE)
end

function distance(p0::pos, p1::pos)
    return sqrt(distance2(p0,p1))
end


function project(p0::pos, distance, bearing)
  dN =  distance * cos(bearing)
  dE =  distance * sin(bearing)
  return pos(p0.n + dN, p0.e + dE, p0.d);
end


function +(p0::pos, p1::pos)
    return pos(p0.n + p1.n , p0.e + p1.e, p0.d + p1.d)
end

#################################################
#Parametrize things with a bearing and a length,
#then we will populate the positions as a tree
#################################################

RefLength = 2000/3;

psi_L = Dict([:T], [(0.5, 180)])
psi_L[:R] = (1, 90)
psi_L[:U1] = (5, 90);
                psi_L[:LX1] = (2, 0);
psi_L[:LD1] = (6, 270); psi_L[:LD2] = (5, 270)
                psi_L[:LB1] = (2, 180);
psi_L[:F1] = (5, 90); #Back on runway


psi_L[:F0] = (3, 90);
psi_L[:GO] = (1, 45);


psi_L[:U2] = (3, 90);
          psi_L[:LX2] = (2, 0)
psi_L[:LD0] = (3, 270); psi_L[:LD3] = (3, 270);
          psi_L[:LB2] = (2, 180);

psi_L[:LDep] = (10, 45);
psi_L[:LArr] = (10, 270);



posNE = Dict([(:R,"S")], [pos(0,0,-300)])
i = 0;
while(true)
  i+=1;
  newleg = false
  for s in allstates
    s_start = (s, "S");   s_end = (s, "E")
    if haskey(posNE, s_start) && !haskey(posNE, s_end)
      posNE[s_end] = project(posNE[s_start],
                             psi_L[s][1]*RefLength,
                             deg2rad(psi_L[s][2]))
      for snext in NextStates[s]
        sn_start = (snext, "S")
        if !haskey(posNE, sn_start) && haskey(psi_L, snext)
          posNE[sn_start] = deepcopy(posNE[s_end])
          newleg = true
        end
      end
    end
  end
  if newleg == false
    break
  end
end
for k in keys(posNE)
  if k[1] in [:R, :T] || k == (:F1, "E")
    posNE[k].d = 0
  elseif k[1] in [:LB1, :F0] && k[2] == "E"
    posNE[k].d = -100
  end
end
posNE[(:F1, "S")].d = -100;


allstates_string = [string(a) for a in allstates]
for astr in allstates_string
  if(astr[1] == 'R')
    astr_l = replace(astr, 'R', 'L')

    for d in {"S", "E"}
      a = (symbol(astr), d)
      b = (symbol(astr_l), d)
      if b in keys(posNE)
        posNE[a] = deepcopy(posNE[b])
        posNE[a].n *= -1;
      end
    end
  end
end


function wpPos(wp::Symbol)
  return deepcopy(posNE[(wp, "S")])
end

#################################################
type airplane
#################################################
  #State:
  airspeed::Float64
  posNED::pos

  psi::Float64
  roll::Float64
  gamma::Float64

  #Keep track of whether we are ready to transition
  #And ready for an ATC command
  readyToTransition::Bool
  readyForATC::Bool

  #VS1 stores the initialization speed of the A/C
  VS1::Float64
  #Noise in the navigation
  navNoise::pos

  #Destination we're heading towards
  navDest::(Symbol, String)
  destNED::pos
  #Command if we've received any?
  atcCommand::Symbol

  #this just keeps track of the history!
  path::Vector{pos}
  sLocHist::Vector{Symbol}


  #Constructors
  function airplane(airspeed, s)
    navOrig = (s, "S")
    navDest = (s, "E")
    p0 = deepcopy(posNE[navOrig])
    p1 = deepcopy(posNE[navDest])
    psi = bearing(p0, p1)
    destNED = p1

    new(airspeed, p0, psi, 0, 0, false, false,
        airspeed, pos(0,0,0),
        navDest, destNED,
        :∅, [deepcopy(p0)], [s])

  end
  airplane(airspeed) = airplane(airspeed, :R)

end

#Find out where this airplane is headed
#(accounting for noise)
function destination(a::airplane)
  return a.destNED #posNE[a.navDest] + a.navNoise
end


#"Rigid" body dynamics of 3DOF sim
#################################################
function move!(ac::airplane, dt::Float64, savepath::Bool = true)
#################################################
  #Slow down if we are taxiing
  if ac.navDest[1] == :T
    ac.airspeed = taxiSpeed
  #Accelerate on the runway to takeoff
  elseif ac.navDest[1] == :R && ac.navDest[2] == "E"
    ac.airspeed = min(ac.airspeed + 1, ac.VS1)
  end

  #Euler step for Position
  dN =  ac.airspeed * cos(ac.psi) * dt
  dE =  ac.airspeed * sin(ac.psi) * dt
  dD = -ac.airspeed * sin(ac.gamma) * dt
  ac.posNED.n += dN;
  ac.posNED.e += dE;
  ac.posNED.d += dD;

  #Coordinated turn
  psidot = 9.81 * tan(ac.roll) / ac.airspeed
  ac.psi = unwrap(ac.psi + psidot * dt);

  if(savepath)
    push!(ac.path,deepcopy(ac.posNED))
    push!(ac.sLocHist, ac.navDest[1])
  end
  #push!(ac.psiHist,ac.psi)
  #push!(ac.rollHist,ac.roll)

end

#################################################
function aviate!(ac::airplane, altitude_desired::Float64, heading_desired::Float64)
#################################################
  #Climb towards the desired altitude
  const r2d = 180./pi
  altitude = -ac.posNED.d

  kp_g =  10. /100.
  if(altitude_desired < 1. && altitude < 2.)
    kp_g *= 10.
  end
  ac.gamma = deg2rad(clip( kp_g * (altitude_desired - altitude),
                          -5., 5.))

  #Roll controller

  #Profiling, putting each on its own line
  heading_error_deg = unwrap(heading_desired - ac.psi)*r2d #+ randn(rng)*4.;

  const kp_roll = 3.0 #clip(randn(rng) + 1, 0.5, 1.5) * 2
  ac.roll = clip( heading_error_deg * kp_roll, -45., 45.)/r2d


  #Special case on the ground,
  #just point in the heading we want directly!
  if(ac.airspeed <= taxiSpeed)
    ac.psi = heading_desired
    ac.roll = 0.
  end
end


#################################################
function navigate!(ac::airplane)
#################################################
  p0 = ac.posNED
  p1 = destination(ac)

  #We will signal we are ready for transition
  #if the distance to the target is below a threshold
  #If the target is @ "E", we are also ready to transition
  ac.readyToTransition = (distance(p1, ac.posNED) < transThresh)
  ac.readyForATC = ac.readyToTransition && ac.navDest[2] == "E"

  #When tracking the runway, do more of a fake x-track like
  if(ac.navDest[1] == :F1 && ac.navDest[2] == "E" && abs(p0.e) < 2000 )
    frac = clip(abs(p0.e)/1500, 0., 1.)
    p1 = pos(p1.n, (p0.e + 100) * (1-frac) + p1.e * frac , p1.d)
  end

  #After navigating, we should aviate
  aviate!(ac, -p1.d, bearing(p0, p1))

end



#################################################
function transition(ac::airplane)
#################################################

  #Where were we heading
  s = ac.navDest[1]

  #If it's towards the start of a leg
  #We will transition to the End point
  if(ac.navDest[2] == "S")
    ac.navDest = (s, "E")

    #This is where we inject some noise to make
    #Things more 'realistic'. Except, No Noise on the runway!
    if ac.navDest == (:F1, "E") || s == :R || s == :T
      ac.navNoise = pos(0,0,0)
    else
      ned = randn(rng,3) .* Float64[maxNoise, maxNoise, maxNoiseAlt]
      ac.navNoise = pos(ned...)
    end
  else
    #If we arrived to the end of a leg, we need to decide where to go next
    #We do that based on any atcCommand that we have received
    a = ac.atcCommand
    sn = randomChoice(s, a != :∅, a)
    ac.atcCommand = :∅


    #Special handling of the departure state
    #If we are departed and we are staying departed
    #Don't go back to the start point, instead linger
    #around the end point for 30 seconds
    d = "S"
    if (s == sn && (s == :LDep || s == :RDep))
      d = "E"
      ac.navNoise = project(pos(0,0,0),
                            30 * ac.airspeed + transThresh,
                            rand(rng)*(2*pi));
    #Also for the go around state, we should head straight
    #to the go around state
    elseif (s == :GO)
      d = "E"
    end
    ac.navDest = (sn, d)

  end

  ac.destNED = posNE[ac.navDest] + ac.navNoise


end

#################################################
function flyPattern!(ac::airplane)
#################################################
  #Check if we are ready to transition based
  #on the last navigation step. If so act accordingly
  if(ac.readyToTransition || ac.atcCommand == :GO)
    transition(ac)
  end

  #Afterward go and navigate
  #(moving will be done later)
  navigate!(ac)
end


#################################################
function runAutoATC(acList::Vector{airplane}, runATC::Symbol, policyFun)
#################################################
  act = g_noaction
  if runATC == :MDP
    act = policyFun([ac.navDest[1] for ac in acList])
#   else
#     for i in 1:4
#       ac = acList[i]
#       if(ac.navDest[1] == :T && ac.navDest[2] == "E")
#         allok = true
#         for j in [1:(i-1) , (i+1):4]
#           ac2 = acList[j]
#           if(ac2.navDest[1] == :R || (ac2.navDest[1] == :F1 && ac2.navDest[2] == "E"))
#             allok = false
#             break
#           end
#         end
#         if allok
#           act = (i, :R)
#         end
#       end
#     end
  end
  return act
end



function isSafe(s::Symbol, dest::ASCIIString)
  #safeStates = Symbol[:T, :LDep , :RDep, :LArr, :RArr]
  return s == :T || s == :LDep || s == :RDep || s == :LArr || s == :RArr || (s == :R && dest == "S")
  #return (s in safeStates)
#   for s2 in safeStates
#     if s == s2
#       return true
#     end
#   end
#   return false
end
#################################################
function getDmin(acList::Vector{airplane})
#################################################
#Compute the minimum distance to a given aircraft
  dmin = Inf
  idmin = Int64[0, 0]
  for idx in 1:(length(acList)-1)
    ac = acList[idx]
    #All is good if we are in a safe state
    if !isSafe(ac.navDest[1], ac.navDest[2])
      #Otherwise iterate over the other aircraft
      #And find the closest one
      for i in (idx+1):length(acList)
        ac2 = acList[i]
        if !isSafe(ac2.navDest[1], ac2.navDest[2])
          dmin_new = distance2(ac.posNED, ac2.posNED)
          if(abs(ac.posNED.d - ac2.posNED.d) < 30 && dmin_new < dmin)
            dmin = dmin_new
            idmin = Int64[idx, i]
          end
        end
      end
    end
  end

  #If we don't find anything, make it NaN for
  #making it straightfoward
  if dmin == Inf
      dmin = NaN
  else
      dmin = sqrt(dmin) - 150.
  end
  return (dmin, idmin)
end

#################################################
function simulate!(acList::Vector{airplane}, Tend, stopEarly = false, runATC::Symbol = :MDP, savepath = true,  policyFun = PiStarFunSlow)
#################################################
#Running simulation,

  #Total time range
  trange =  0:simdt:Tend

  smartATCtiming = (runATC == :MDP || runATC == :None)
  if runATC == :MDP2
    runATC = :MDP
  end


  stopsim = false;
  idmin = Int64[0,0]
  alertCount = 0

  flightTime = 0.
  tidx = 0
  for t in trange
    tidx += 1

    #Find out if any of the aircraft in the pattern
    #is ready for a command. This could be done more
    #concisely but list comprehensions seem to slow
    #things down!

    #Find out if anyone is busy
    noPendingCommand = true
    for idx in 1:length(acList)
      noPendingCommand = noPendingCommand && (acList[idx].atcCommand == :∅)
    end
    noPendingCommand = true

    readyForCommand = false
    if(noPendingCommand) #Someone is busy, don't send a new command
        if(smartATCtiming)
          #Find out if any aircraft is ready for a command
          for idx in 1:length(acList)
            readyForCommand = readyForCommand || acList[idx].readyForATC
          end
        else
          #Do it based on clock
          readyForCommand = (tidx % 40 == 0)
        end
    end

    #If any aircraft is about to transition,
    #see if there's an ATC command that should be issued!
    if(readyForCommand)
      act = runAutoATC(acList, runATC, policyFun)
      #If we have an action to issue, pass it along
      if act != g_noaction && acList[act[1]].atcCommand == :∅
        acList[act[1]].atcCommand = act[2]
        alertCount += 1
      end
    end

    #Fly pattern logic for all aircraft
    for idx in 1:length(acList)
      flyPattern!(acList[idx])
      move!(acList[idx], simdt, savepath)
    end


    #Compute the distance to all other boogies
    (dmin, idmin_local) = getDmin(acList)
    if(stopEarly && dmin <= 0)
      idmin = idmin_local
      break;
    end




    for idx in 1:length(acList)
      if acList[idx].navDest[1] != :T
        flightTime += simdt
      end
    end

#     allInTaxi = true
#     for idx in 1:length(acList)
#       allInTaxi = allInTaxi && acList[idx].navDest[1] == :T && acList[idx].navDest[2] == "S"
#     end
#     if(allInTaxi)
#         tidx = -1;
#         break
#     end

  end

  tmax = 1e6
  if(tidx != -1)
    tmax = trange[tidx]
  end
  return (idmin, tmax, alertCount, flightTime)
end



#################################################
function dummySim(ac::airplane, Dt)
#################################################
    a = true; b = true;
    for t in 0:simdt:Dt
#         if(t < Dt/4)
#             #ac.roll = 0
#             aviate!(ac, 300, ac.psi)
#         elseif(t < Dt/2)
#             #ac.roll = pi/10;
#             aviate!(ac, -ac.posNED.d, deg2rad(90))
#         else
#             #ac.roll = -pi/10;
#             aviate!(ac, -ac.posNED.d, deg2rad(190))
#         end
#         navigate!(ac)

        flyPattern!(ac)

        move!(ac, simdt)
    end
end
