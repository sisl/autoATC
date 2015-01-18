import Base.print
import Base.show
import Base.string

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

simdt=0.25
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



function +(p0::pos, p1::pos)
    return pos(p0.n + p1.n , p0.e + p1.e, p0.d + p1.d)
end

#################################################
#Parametrize things with a bearing and a length, then we
#will populate the positions as a tree
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


function project(p0::pos, distance, bearing)
  dN =  distance * cos(bearing)
  dE =  distance * sin(bearing)
  return pos(p0.n + dN, p0.e + dE, p0.d);
end


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


  VS1::Float32
  navNoise::pos

  navDest::(Symbol, String)
  atcCommand::Symbol

  #this just keeps track of the history!
  path::Vector{pos}



  function airplane(airspeed, s)
    navOrig = (s, "S")
    navDest = (s, "E")
    p0 = deepcopy(posNE[navOrig])
    p1 = deepcopy(posNE[navDest])
    psi = bearing(p0, p1)

    new(airspeed, p0, psi, 0, 0,
        airspeed, pos(0,0,0),
        navDest, :∅, [deepcopy(p0)])

  end
  airplane(airspeed) = airplane(airspeed, :R)

end
function string(r::airplane)
    return string("pos = ", string(r.posNED))
end
print(io::IO, x::airplane) = print(io, string(x))
show(io::IO, x::airplane) = print(io, x)


function destination(a::airplane)
  return posNE[a.navDest] + a.navNoise
end


taxiSpeed = 8
function move(ac::airplane, dt)
  if ac.navDest[1] == :T
    ac.airspeed = taxiSpeed
  elseif ac.navDest == (:R, "E")
    ac.airspeed = min(ac.airspeed + 1, ac.VS1)
  end
  dN =  ac.airspeed * cos(ac.psi) * dt
  dE =  ac.airspeed * sin(ac.psi) * dt
  dD = -ac.airspeed * sin(ac.gamma) * dt


  ac.posNED.n += dN;
  ac.posNED.e += dE;
  ac.posNED.d += dD;

  psidot = 9.81 * tan(ac.roll) / ac.airspeed
  ac.psi = unwrap(ac.psi + psidot * dt);

  push!(ac.path,deepcopy(ac.posNED))
  #push!(ac.psiHist,ac.psi)
  #push!(ac.rollHist,ac.roll)

end


function aviate(ac::airplane, altitude_desired, heading_desired)
    kp =  100 * pi/180
    ac.gamma = deg2rad(clip( kp * (altitude_desired - -ac.posNED.d),
                            -8, 8))

    heading_desired = unwrap(heading_desired)
    heading_error_deg = rad2deg(unwrap(heading_desired - ac.psi)) +
                        randn()*4;

    kp = clip(randn() + 1, 0.5, 1.5) * 2
    ac.roll = deg2rad( clip( heading_error_deg * kp, -45, 45))


    #on the groun, just point in the heading we want directly!
    if(ac.airspeed <= taxiSpeed)
      ac.roll = 0
      ac.psi = heading_desired
    end

end


function navigate(ac::airplane)
        p0 = ac.posNED
        p1 = destination(ac)


        if(ac.navDest == (:F1, "E"))
            p1.e = p0.e + 100
        end

        aviate(ac, -p1.d, bearing(p0, p1))

end



transThresh = 150
maxNoise = 200 #FIXME: Make this configurable per airplane?
function transition(ac::airplane)
  #Assume we've arrived to our destination, move to next waypoint
  #For now deterministic...

  s = ac.navDest[1]
  if(ac.navDest[2] == "S")
    ac.navDest = (s, "E")


    #No Noise on the runway!
    if ac.navDest == (:F1, "E") || s == :R || s == :T
      ac.navNoise = pos(0,0,0)
    else
      ned = [clip(x * maxNoise, -maxNoise, maxNoise) for x in randn(3)]
      ned[3] = 0#clip(ned[3], -maxNoise / 10, maxNoise/10)
      ac.navNoise = pos(ned...)
    end
  else
    sn = randomChoice(s, false, :∅) #Fixme! Use actual atc command!

    d = "S"
    #Special handling of the departure state
    if (s == sn && (s == :LDep || s == :RDep))
      d = "E"
    elseif (s == :GO)
      d = "E"
    end
    ac.navDest = (sn, d)

  end



end

function flyPattern(ac::airplane)
      #Check if we have arrived to our destination
      #If so, move to next waypoint!

       if(distance(destination(ac), ac.posNED) < transThresh)
         transition(ac)
       end

      navigate(ac)
end

function dummySim(ac::airplane, Dt)
    a = true; b = true;
    for t in 0:simdt:Dt
#         if(t < Dt/4)
#             #ac.roll = 0
#             aviate(ac, 300, ac.psi)
#         elseif(t < Dt/2)
#             #ac.roll = pi/10;
#             aviate(ac, -ac.posNED.d, deg2rad(90))
#         else
#             #ac.roll = -pi/10;
#             aviate(ac, -ac.posNED.d, deg2rad(190))
#         end
#         navigate(ac)

        flyPattern(ac)

        move(ac, simdt)
    end
end

function simulate(acList::Vector{airplane}, Dt)
  for t in 0:simdt:Dt
    for ac in acList
      flyPattern(ac)
      move(ac, simdt)
    end
  end
end
