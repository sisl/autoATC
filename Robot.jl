module Robot

export vec2, action, robot, move

using Distributions
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

function Base.isequal(v1::vec2, v2::vec2)
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
  #this just keeps track of the history!
  path::Vector{vec2}
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

#################################################
function move(r::robot, a::action)
#################################################
    #Move based on 'stochastic' process
    γ = a.issuer == auto ? 0.1 : 0.2
    d = Bernoulli(γ)
    α = rand(d); β = rand(d)
    #If we are not healthy, fall towards the ground!
    w_z = r.health == healthy ? 0. : (a.issuer == auto ? -2. : -1.)
    r.position.x += a.velCmd.x * (1-α)
    r.position.z += a.velCmd.z * (1-β) + w_z
    r.position.z = max(r.position.z , 0) #don't go below the ground!


    popTask = length(r.tasks) > 0 && r.tasks[1] == r.position
    if(popTask)
        r.tasks = r.tasks[2:end]
    end

    push!(r.path, vec2(r.position.x, r.position.z))
    return r
end


end