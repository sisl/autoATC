module RewardFun

using auxFuns
using pattern

export Reward

#It is assumed that some other module (pattern) is providing the following:

# g_allStates: an array of symbols representing all possible substates that we can be in
# g_sn: a dictionary of type (Symbol => Int64) which goes from substate to index
# legalActions(S): function which tells us which actions are allowed!
# g_nullAct: the noaction case


###########################################
###########################################
###########################################
###########################################

###########################################
#Defining Reward functions
###########################################
#Each collision costs 1000.
#And each aircraft just sitting on the taxi also incurs cost
const collisionCost = -1000.0f0
const taxiCost = -10.0f0


#############################################
function Reward(SX::Union(XType, SType), a::ActType, β::Float32; 
                timeHorizon::Float32 = 1f0, 
                E::Vector{Float32} = zeros(Float32, pattern.g_nVehicles)
                )
#############################################
    R = 0.0f0

    #Actions have a cost
    if (!pattern.isNullAct(a))
        R += β * collisionCost;
    end
    
    #Check for collisions
    R += NcolNtaxi(SX, collisionCost, taxiCost, timeHorizon, E)
    
    return R;
end


###########################################
function NcolNtaxi(SX::Union(XType, SType), 
                  collisionCost::Float32, 
                  taxiCost::Float32,
                  timeHorizon::Float32,
                  E::Vector{Float32}
                  )
###########################################
  Nc = 0f0
  Nt = 0
  for i in 1:length(SX)
    isTaxi = false
    isSafe = false
    if typeof(SX) == XType
        isTaxi = isin(SX[i], pattern.xTaxi)
        isSafe = isin(SX[i], pattern.xSafe)
    else
        isTaxi = isin(SX[i], pattern.sTaxi)
        isSafe = isin(SX[i], pattern.sSafe)
    end
    #Count vehicles in Taxi
    if(isTaxi)
        Nt += 1
    end
    #Count vehicles in collision
    if isSafe #Ignore vehicles that are safe
      continue
    end
    for j in (i+1):length(SX)
      if SX[j] == SX[i]
        dt = abs(E[i] - E[j])
        if dt <= timeHorizon
            Nc += 1f0
        else
            Nc += max(0f0, 1f0 - (dt-timeHorizon)/timeHorizon/3f0)
        end
      end
    end
  end
  return Nc*collisionCost +  Nt*taxiCost
end

end