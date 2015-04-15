module airplaneType

using posType
using pattern

export airplane, move!, flyPattern!, getDmin!

#Parameters
type simParsType
    transThresh::Float64
    maxNoise::Float64 #FIXME: Make this configurable per airplane? ~500 ft
    maxNoiseAlt::Float64 #FIXME: Make this configurable per airplane? ~60 ft
    taxiSpeed::Float64 #default taxi speed
    
    function simParsType()
        new(50., #transThresh
           200., #maxNoise ~500ft
           20.,  #maxNoise altitude ~60ft
           5., #taxi speed
           )
    end
end
const simPars=simParsType()


const DebugOn = false
#################################################
#Parametrize things with a bearing and a length,
#then we will populate the positions as a tree
#################################################
const RefLength = 2000/3;

rhotheta = (dn, de) -> (sqrt(de^2+dn^2), rad2deg(atan2(de,dn)))

psi_L = Dict([:T], [(0.5, 180.)])
psi_L[:R] = (1, 90)
psi_L[:U1] = (5, 90);
                psi_L[:LX1] = (2, 0);
psi_L[:LD1] = (6, 270); psi_L[:LD2] = (5, 270)
                psi_L[:LB1] = (2, 180);
psi_L[:F1] = (4, 90); #Back on runway


psi_L[:F0] = (3, 90);
psi_L[:GO] = rhotheta(0.3, 2) #(1.5, 40);


psi_L[:U2] = (3, 90);
          psi_L[:LX2] = (2, 0)
psi_L[:LD0] = (3, 270); psi_L[:LD3] = (3, 270);
          psi_L[:LB2] = (2, 180);

psi_L[:LDep] = (10, 45);
psi_L[:LArr] = rhotheta(-3, -16) #3 south, 16 west



posNE = Dict([(:R,"S")], [pos(0,0,-300)])
while(true)
  newleg = false
  for s in g_allstates
    s_start = (s, "S");   s_end = (s, "E")
    if haskey(posNE, s_start) && !haskey(posNE, s_end)
      posNE[s_end] = project(posNE[s_start],
                             psi_L[s][1]*RefLength,
                             deg2rad(psi_L[s][2]))
      for snext in NextStates[s]
        snext = phaseFree(snext)
        sn_start = (snext, "S")
        if !haskey(posNE, sn_start) && haskey(psi_L, snext)
          posNE[sn_start] = copy(posNE[s_end])
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


for astr in g_allstates_string
  if(astr[1] == 'R')
    astr_l = replace(astr, 'R', 'L')

    for d in {"S", "E"}
      a = (symbol(astr), d)
      b = (symbol(astr_l), d)
      if b in keys(posNE)
        posNE[a] = copy(posNE[b])
        posNE[a].n *= -1;
      end
    end
  end
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
  destHist::Vector{pos}
  noiseHist::Vector{pos}
  #Leg distance
  legDist::Float64
  navPhase::Int64


  #Constructors
  function airplane(airspeed, s, frac)
    s = phaseFree(s)
    navOrig = (s, "S")
    navDest = (s, "E")
    p0 = copy(posNE[navOrig])
    p1 = copy(posNE[navDest])
    psi = bearing(p0, p1)
    destNED = p1
    
    p0.n = (1-frac) * p0.n + frac* p1.n
    p0.e = (1-frac) * p0.e + frac* p1.e
    p0.d = (1-frac) * p0.d + frac* p1.d

    legDist = distance(p0,p1)
    new(airspeed, p0, psi, 0, 0, false, false,
        airspeed, pos(0,0,0),
        navDest, destNED,
        :∅,
        [copy(p0)], [s], [copy(p1)], [pos()],
        legDist, 1)

  end
  airplane(airspeed, s) = airplane(airspeed, s, 0.)
  airplane(airspeed) = airplane(airspeed, :R, 0.)

end


#"Rigid" body dynamics of 3DOF sim
#################################################
function move!(ac::airplane, dt::Float64, savepath::Bool = true)
#################################################
  #Slow down if we are taxiing
  if ac.navDest[1] == :T
    ac.airspeed = simPars.taxiSpeed
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
    push!(ac.path,copy(ac.posNED))
    push!(ac.sLocHist, ac.navDest[1])
  end

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
  if(ac.airspeed <= simPars.taxiSpeed)
    ac.psi = heading_desired
    ac.roll = 0.
  end
end


#################################################
function navigate!(ac::airplane)
#################################################
  p0 = ac.posNED
  p1 = ac.destNED

  #We will signal we are ready for transition
  #if the distance to the target is below a threshold
  #If the target is @ "E", we are also ready to transition
  d = distance(p1, p0)
  ac.readyToTransition = (d < simPars.transThresh)
  ac.readyForATC = ac.readyToTransition && ac.navDest[2] == "E"
  
  ac.navPhase = min(max(int(ceil((1 - d/ ac.legDist) * nPhases)), 1) , nPhases)  
  
  #When tracking the runway, do more of a fake x-track like
  if(ac.navDest[1] == :F1 && ac.navDest[2] == "E" && abs(p0.e) < 2000 )
    frac = clip(abs(p0.e)/1500, 0., 1.)
    p1 = pos(p1.n, (p0.e + 100) * (1-frac) + p1.e * frac , p1.d)
  end

  #After navigating, we should aviate
  aviate!(ac, -p1.d, bearing(p0, p1))

end



#################################################
function transition!(ac::airplane)
#################################################

  #Where were we heading
  s = ac.navDest[1]

  #If it's towards the start of a leg
  #We will transition to the End point
  if(ac.navDest[2] == "S")
    ac.navDest = (s, "E")

    #This is where we inject some noise to make
    #Things more 'realistic'. Except, No Noise on the runway!
    if ac.navDest == (:F1, "E") || s == :R
      ac.navNoise = pos(0,0,0)
    else
      ned = randn(pattern.rng,3) .* Float64[simPars.maxNoise, simPars.maxNoise, simPars.maxNoiseAlt]
      
      #Departure states get a bit more noise
      if s == :LDep || s == :RDep
        ned[1] *= 5
        ned[2] *= 5
      #Taxi state gets less noise
      elseif s == :T
        ned[1] *= 0.1
        ned[2] *= 0.1
      end
      
      ac.navNoise = pos(ned)
    end
  else
    #If we arrived to the end of a leg, we need to decide where to go next
    #We do that based on any atcCommand that we have received
    a = ac.atcCommand
    sn = randomChoice(s, a != :∅, a)
    #Get rid of the phase information
    sn = phaseFree(sn)
    ac.atcCommand = :∅


    #Special handling of the departure state
    #If we are departed and we are staying departed
    #Don't go back to the start point, instead linger
    #around the end point for 30 seconds
    d = "S"
    if (s == sn && (s == :LDep || s == :RDep))
      d = "E"
      ac.navNoise = project(pos(0,0,0),
                            30 * ac.airspeed + simPars.transThresh,
                            rand(pattern.rng)*(2*pi));
    #Also for the go around state, we should head straight
    #to the end of the leg!
    elseif (s == :GO)
      d = "E"
    end
    ac.navDest = (sn, d)

  end

  ac.destNED = posNE[ac.navDest] + ac.navNoise
  
  if(ac.navDest == (:GO, "E"))
    #Go around state is either to left or to right
    ac.destNED.n *= rand(pattern.rng) > 0.5 ? -1. : 1.
  end

  #Compute total distance to be travelled. This will be used
  #to guess what phase the aircraft are in. 
  ac.legDist = distance(ac.destNED, ac.posNED)

  #TODO: Consider putting navPhase as part of navDest?
  #also reset the navPhase
  ac.navPhase = 1 
  
  #Don't waste navigating towards the start point if it's right next
  #to where we are going!
  if(ac.legDist <= simPars.transThresh)
    transition!(ac)
    return;
  end
  
  if(DebugOn)
    push!(ac.destHist, copy(ac.destNED))
    push!(ac.noiseHist, copy(ac.navNoise))
  end
end

#################################################
function flyPattern!(ac::airplane, simdt::Float64, savepath::Bool=true)
#################################################
  #Check if we are ready to transition based
  #on the last navigation step. If so act accordingly
  if(ac.readyToTransition) # || phaseFree(ac.atcCommand) == :GO)
    transition!(ac)
  end

  #Afterward go and navigate
  #(moving will be done later)
  navigate!(ac)
  
  move!(ac, simdt, savepath)

end



#################################################
function isSafe(s::Symbol, dest::ASCIIString)
#################################################
  return s == :T || s == :LDep || s == :RDep || s == :LArr || s == :RArr || (s == :R && dest == "S")
end
#################################################
function getDmin!(idmin, acList::Vector{airplane})
#################################################
#Compute the minimum distance to a given aircraft
  dmin = Inf
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
            idmin[1] = idx; idmin[2] = i
          end
        end
      end
    end
  end

  #If we don't find anything, make it NaN for
  #making it straightfoward
  if dmin != Inf
      dmin = sqrt(dmin) - 150.
  end
  return dmin
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

end #end of module

