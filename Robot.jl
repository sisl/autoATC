#module Robot

#export vec2, action, robot, move, simulate, z

# macro quickuse(mod)
#   if isdefined(current_module(),:LastMain) && isdefined(LastMain,:$mod)
#     ex = :(using LastMain.$mod)
#   else
#     ex = :(using $mod)
#   end
# end
# @quickuse Distributions
# @quickuse PyPlot

if isdefined(current_module(),:LastMain) && isdefined(LastMain,:PyPlot)
  using LastMain.PyPlot
else
  using PyPlot
end

if isdefined(current_module(),:LastMain) && isdefined(LastMain,:Distributions)
  using LastMain.Distributions
else
  using Distributions
end



import Base.isequal
import Base.print
import Base.show
import Base.string
#################################################
immutable pilotType
#################################################
  val::Int8
end
auto = pilotType(0)
human = pilotType(1)
#################################################
immutable attAction
#################################################
  val::Int8
end
const set0 = attAction(0)
const set1 = attAction(1)
const set2 = attAction(2)
const up = attAction(3)
const dn = attAction(4)
const nop = attAction(5)
#################################################
immutable attLvl
#################################################
  val::Int8
end
const attLvl0 = attLvl(0)
const attLvl1 = attLvl(1)
const attLvl2 = attLvl(2)
#################################################
immutable healthType
#################################################
    val::Int8
end
const healthy = healthType(0)
const notHealthy = healthType(1)

#################################################
type vec2
#################################################
  x::Float64 #X position
  z::Float64 #Z position
end
print(io::IO, x::vec2) = print(io, string(x))
show(io::IO, x::vec2) = print(io, x)

function ==(v1::vec2, v2::vec2)
    v1.x == v2.x && v1.z == v2.z
end

function string(v::vec2)
    return string("(",string(v.x),",",string(v.z),")")
end

#################################################
type robot
#################################################
  #State:
  position::vec2 #robot position
  health::healthType    #robot health.
  attention::attLvl #operator's attention {0,1,2}
  tasks::Vector{vec2} # list of tasks
  #Some parameters
  Score::Float64
  crashed::Bool
  γ::Vector{Float64} #Bernoulli parametrizatin of accuracy {robot, human}
  #this just keeps track of the history!
  path::Vector{vec2}

  robot(p,h,a,t,γ) = new(p,h,a,t,0., false, γ, [vec2(p.x, p.z)])
  robot(p,h,a,t) = robot(p,h,a,t,[0.2, 0.1])
end
function string(r::robot)
    return string("pos = ", string(r.position))
end
print(io::IO, x::robot) = print(io, string(x))
show(io::IO, x::robot) = print(io, x)

#################################################
type action
#################################################
    issuer::pilotType
    velCmd::vec2
    attCmd::attAction
end

function clip(x::Number, min::Number, max::Number)
  if(x < min)
    return min
  elseif(x > max)
    return max
  else
    return x
  end
end

#################################################
function move(r::robot, a::action)
#################################################
  #Move based on 'stochastic' process
  if (! r.crashed  )
      γ = a.issuer == auto ? r.γ[1] : r.γ[2]
      d = Bernoulli(γ)
      α = rand(d); β = rand(d)
      #If we are not healthy, fall towards the ground!
      w_z = r.health == healthy ? 0. : (a.issuer == auto ? -2. : -1.)
      r.position.x += clip(a.velCmd.x,-1.,1.) * (1-α)
      r.position.z += clip(a.velCmd.z,-1.,1.) * (1-β) + w_z
      if(r.position.z <= 0 && r.position.x != 0)
        r.crashed = true
        r.Score -= 1e3
      else
        r.Score -= (abs(a.velCmd.x) + abs(a.velCmd.z))
      end
      r.position.z = max(r.position.z, 0)
      push!(r.path, vec2(r.position.x, r.position.z))
  end

  popTask = length(r.tasks) > 0 && r.tasks[1] == r.position
  if(popTask)
    r.tasks = r.tasks[2:end]
    r.Score += 100
  end

  return r
end

#################################################
function simulate(r::robot, policy, randFailure = false, N=100)
#################################################
  while length(r.path) < N && !r.crashed
    #Move according to policy
    move(r,policy(r))

    #Put some spice in life
    if(randFailure)
      if (rand() > 0.99)
        r.health = notHealthy
      end
    elseif(r.position.z > 6)
        r.health = notHealthy
    end

    #No more tasks and we are back home!
    if((length(r.tasks) == 0 || r.health == notHealthy) && r.position == vec2(0,0))
      break
    end
end

end
#end
