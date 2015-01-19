# import Base.print
# import Base.show
# import Base.string

function unwrap(a)
    if (a > pi)
        return mod(a + pi, 2*pi) - pi
    elseif (a < -pi)
        return -unwrap(-a)
    else
        return a;
    end
end

function clip(r, min, max)
  if r < min
    return min
  elseif r > max
    return max
  else
    return r
  end
end


#Parameters
simdt=0.25
transThresh = 150
maxNoise = 200 #FIXME: Make this configurable per airplane?
taxiSpeed = 8 #default taxi speed

#################################################
type pos
#################################################
  n::Float32
  e::Float32
  d::Float32
end

function bearing(p0::pos, p1::pos)
    dN = p1.n - p0.n;
    dE = p1.e - p0.e;
    return unwrap(pi/2-atan2(dN, dE))
end


function distance(p0::pos, p1::pos)
    dN = p0.n - p1.n;
    dE = p0.e - p1.e;
    return sqrt(dN*dN + dE*dE)
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
  airspeed::Float32
  posNED::pos

  psi::Float32
  roll::Float32
  gamma::Float32

  #Keep track of whether we are ready to transition
  #And ready for an ATC command
  readyToTransition::Bool
  readyForATC::Bool

  #VS1 stores the initialization speed of the A/C
  VS1::Float32
  #Noise in the navigation
  navNoise::pos

  #Destination we're heading towards
  navDest::(Symbol, String)
  #Command if we've received any?
  atcCommand::Symbol

  #this just keeps track of the history!
  path::Vector{pos}


  #Constructors
  function airplane(airspeed, s)
    navOrig = (s, "S")
    navDest = (s, "E")
    p0 = deepcopy(posNE[navOrig])
    p1 = deepcopy(posNE[navDest])
    psi = bearing(p0, p1)

    new(airspeed, p0, psi, 0, 0, false, false,
        airspeed, pos(0,0,0),
        navDest, :∅, [deepcopy(p0)])

  end
  airplane(airspeed) = airplane(airspeed, :R)

end

#Find out where this airplane is headed
#(accounting for noise)
function destination(a::airplane)
  return posNE[a.navDest] + a.navNoise
end


#"Rigid" body dynamics of 3DOF sim
#################################################
function move!(ac::airplane, dt, savepath = true)
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
  end
  #push!(ac.psiHist,ac.psi)
  #push!(ac.rollHist,ac.roll)

end

#################################################
function aviate!(ac::airplane, altitude_desired, heading_desired)
#################################################
  #Climb towards the desired altitude
  kp =  100 * pi/180
  ac.gamma = deg2rad(clip( kp * (altitude_desired - -ac.posNED.d),
                          -8, 8))

  #Roll controller
  heading_desired = unwrap(heading_desired)
  heading_error_deg = rad2deg(unwrap(heading_desired - ac.psi)) +
    randn()*4;

  kp = clip(randn() + 1, 0.5, 1.5) * 2
  ac.roll = deg2rad( clip( heading_error_deg * kp, -45, 45))


  #Special case on the ground,
  #just point in the heading we want directly!
  if(ac.airspeed <= taxiSpeed)
    ac.roll = 0
    ac.psi = heading_desired
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
  if(ac.navDest[1] == :F1 && ac.navDest[2] == "E")
    p1.e = p0.e + 100
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
      ned = [clip(x * maxNoise, -maxNoise, maxNoise) for x in randn(3)]
      ned[3] = 0#clip(ned[3], -maxNoise / 10, maxNoise/10)
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
                            rand()*(2*pi));
    #Also for the go around state, we should head straight
    #to the go around state
    elseif (s == :GO)
      d = "E"
    end
    ac.navDest = (sn, d)

  end


end

#################################################
function flyPattern!(ac::airplane)
#################################################
  #Check if we are ready to transition based
  #on the last navigation step. If so act accordingly
  if(ac.readyToTransition)
    transition(ac)
  end

  #Afterward go and navigate
  #(moving will be done later)
  navigate!(ac)
end


#################################################
function runAutoATC(acList::Vector{airplane}, runATC::Symbol)
#################################################
  act = noaction
  if runATC == :MDP
    act = PiStarFunSlow([ac.navDest[1] for ac in acList])
  elseif runATC == :Greedy
    s = [ac.navDest[1] for ac in acList]
    s_i = s2i[s]
    act_i = PiSimple[s_i]
    act = validActions(s)[act_i]
  else
    for i in 1:4
      ac = acList[i]
      if(ac.navDest[1] == :T && ac.navDest[2] == "E")
        act = (i, :R)
        break
      end
    end
  end
  return act
end



safeStates = [:T, :LDep , :RDep, :LArr, :RArr]
#################################################
function getDmin(acList::Vector{airplane}, idx)
#################################################
#Compute the minimum distance to a given aircraft
  dmin = Inf
  #This A/C
  ac = acList[idx]
  #All is good if we are in a safe state
  if !(acList[idx].navDest[1] in safeStates)
    #Otherwise iterate over the other aircraft
    #And find the closest one
    for i in [1:(idx-1) , (idx+1):length(acList)]
      ac2 = acList[i]
      if !(ac2.navDest[1] in safeStates)
        dmin = min(dmin, distance(ac.posNED, ac2.posNED))
      end
    end
  end

  #If we don't find anything, make it NaN for
  #making it straightfoward
  dmin = dmin - 100
  if dmin == Inf
      dmin = NaN
  end
  return dmin
end

#################################################
function simulate!(acList::Vector{airplane}, Tend, stopEarly = false, runATC::Symbol = :MDP, savepath = true)
#################################################
#Running simulation,

  #Total time range
  trange =  0:simdt:Tend

  #Commented out for now
  #RNMAC = zeros(Float32,length(trange), length(acList))

  stopsim = false;

  tidx = 0
  for t in trange
    tidx += 1

    #Find out if any of the aircraft in the pattern
    #is ready for a command
    readyForCommand = any([ac.readyForATC for ac in acList])

    #If any aircraft is about to transition,
    #see if there's an ATC command that should be issued!
    if(readyForCommand)
      act = runAutoATC(acList, runATC)
      #If we have an action to issue, pass it along
      if act != noaction
        acList[act[1]].atcCommand = act[2]
      end
    end

#     Fly pattern logic for all aircraft
    for ac in acList
      flyPattern!(ac)
      move!(ac, simdt, savepath)

      if readyForCommand
         #print(" ");
#           @printf("t = %.2f %s; %s -> %s\n",
#                   2.2, "a","b","c")#ac.posNED, ac.navDest, ac.atcCommand)
      end
    end


    #Compute the distance to all other boogies
    for acidx in 1:4
      dmin = getDmin(acList, acidx)
      if(stopEarly && dmin <= 0)
        stopsim = true
        break;
      end
      #RNMAC[tidx, acidx] = dmin;

    end
    if stopsim
      break
    end
  end

  return (0, trange[tidx])
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
