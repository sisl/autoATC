module CTMDP_kronsolver

using pattern
using RewardFun

using HDF5, JLD


export loadPolicy, savePolicy
export ctmdpPolicy


function combos_with_replacement(list, k)
    n = length(list)
    [[list[c[i]-i+1] for i=1:length(c)] for c in combinations([1:(n+k-1)],k)]
end


#TODO: Find a way to not have to allocate all of this memory?
###########################################
###########################################
#Awesome functions for indexing magic :)
###########################################
#Convention is as follows:
#s is the substate in the pattern of a single vehicle as a symbol
#x is the substate in the pattern of a single vehicle as an index

#S is the state as an array of symbols of all vehicles S = [s1, s2, ...]
#X is the state as an array of Int64 of all vehicles X = [x1, x2, ...]

#CIDX is the compact index of the state (where order does not matter) there are C(n+k-1, k) of them
#LIDX is the long    index of the state (where order does matter)     there are n^k of them
###########################################
###########################################

#Get all possible states, allowing replacement, but order does not matter:
const g_Scomp = combos_with_replacement(g_allstates, g_nVehicles)
const g_Xcomp = XType[S2X(s) for s in g_Scomp]

const g_nScomp = length(g_Scomp); const g_nXcomp = g_nScomp;
const g_nSlong  = g_nNodes^g_nVehicles; const g_nXlong = g_nSlong;

#Going from compact indices to states
function CIDX2S(cidx::Int64)
  return g_Scomp[cidx]
end
function CIDX2X(cidx::Int64)
  return g_Xcomp[cidx]
end

#Going from states to compact indices
const g_X2CIDX_dict = (XType => Int64)[g_Xcomp[cidx] => cidx for cidx in 1:g_nScomp]
function X2CIDX(X::XType)
  return g_X2CIDX_dict[sort(X)]
end
function S2CIDX(S::SType)
  return X2CIDX(S2X(S))
end

#This is when we need to go between compact and long indices
#Note that many LIDX map to the same CIDX
#And that one CIDX maps to many LIDX (but we only return one)
function LIDX2CIDX(lidx::Int64)
  #Actually, this is the one one we really
  #care about making faster. Might need to get rid
  #of the reverse being done?
  X = LIDX2X(lidx);
  return X2CIDX(X);
end
function CIDX2LIDX(cidx::Int64)
  return X2LIDX(g_Xcomp[cidx])
end


#Until we can trust this like crazy, doing some
#random unit tests

using Base.Test


function unitTest(idx, isCompact)
  cidx = 0; lidx = 0
  if(isCompact)
    cidx = idx; lidx = CIDX2LIDX(idx)
  else
    lidx = idx; cidx = LIDX2CIDX(idx)
  end
  Xc = CIDX2X(cidx); Sc = CIDX2S(cidx)
  Xl = LIDX2X(lidx); Sl = LIDX2S(lidx)
  if (!isCompact)
    #In this case, we can't guarantee
    #That the lidx
    Xl = sort(Xl);
    Sl = X2S(Xl)
  end
  @test Xc == Xl == sort(Xl) == sort(Xc)
  @test Xc == S2X(Sl) == S2X(Sc)
  @test Sc == Sl == X2S(Xc) == X2S(Xl)
end
for cidx in rand(1:g_nXcomp, 100)
  unitTest(cidx, true)
end
for lidx in rand(1:g_nXlong, 100)
  unitTest(lidx, false)
end

###########################################
###########################################
###########################################
###########################################



function QVeval(X::XType, action::compActType, Qt_list, V::Vector{Float32},
                ζ::Float32, β::Float32, V_is_compact::Bool, 
                Cb_res, res_u, res_u_rowval, res_u_nzval)

  R = Reward(X, action, β)

  #This is a terminal state... Q(s,a) = R(s,a)
  assert(β < 0.9f0) #We make the assumption that action cost is small relative to collision cost
  if( R <=  collisionCost)
    return R
  end
  

  (idx, act) = action;

  #First, re-order the states so that we are addressing the right instance
  swap!(X, 1, idx)

  #Find the long index
  X_cidx = X_lidx  = X2LIDX(X)
  if(V_is_compact)
    X_cidx = X2CIDX(X)
  end 
  
  
  #get the relevant Q row
  Qt = Qti_ABt(Qt_list[act+1], Qt_list[1], g_nVehicles-1, X_lidx, 
                    Cb_res, res_u, res_u_rowval, res_u_nzval)

  #These are the entries of the transition which are non
  #zero. Note that Qt is a column vector
  Qsparse_indices = Qt.colptr[1] : (Qt.colptr[2]-1)

  #Diagonal element is negative by construction
  #But q(x) is defined as the positive value
  qx = -Qt[X_lidx]

  #q(x) is not supposed to be part of the summation
  #we will later add Qt[si,si]*V[si] so starting with -Qt[si,si] will cancel it out
  qVsum =  qx * V[X_cidx]

  for Qsparse_idx in Qsparse_indices
    Qval = Qt.nzval[Qsparse_idx]
    V_CIDX = V_LIDX = Qt.rowval[Qsparse_idx]
    if(V_is_compact)
        V_CIDX = LIDX2CIDX(V_LIDX)
    end     

    qVsum += Qval * V[V_CIDX]
  end
  
  #Put back X in its original state
  swap!(X, 1, idx)

  return (qVsum ) / (ζ + qx) + R

end


function gaussSeidel!(Qt_list, V::Vector{Float32}, ζ::Float32, β::Float32; maxIters::Int64=100, maxTime::Float64 = Inf)

  Aopt = (compActType)[copy(pattern.g_nullAct) for i in 1:g_nXcomp];
  
  V_is_compact = length(Aopt) == length(V)

  compActs = compActType[copy(pattern.g_nullAct) for i in 1:pattern.g_nMaxActs]

  Xp_indices = collect(permutations(1:g_nVehicles))
  
  n = size(Qt_list[1],1)
  res_u = spzeros(Float32, n,1)
  res_u_rowval = Array(Int64, n)
  res_u_nzval = Array(Float32, n)
  Cb_res = spzeros(Float32, n^(g_nVehicles-1),1); #n^K with K = nVehicles - 1
  
  
  start = time()
  @time for iter in 1:maxIters
    maxVchange = 0.0f0
    X_cidx = 0
    nActsChanged = 0;  
    for X in g_Xcomp
        X_cidx += 1 #X_cidx = X2CIDX(X)
        aopt = pattern.g_nullAct
        Qmax = float32(-Inf)
        #Populate compact actions for this Xtate
        nActs = pattern.validCompActions!(compActs, X)
        for aIdx in 1:nActs
            Qa = QVeval(X, compActs[aIdx], Qt_list, V, ζ, β, 
                        V_is_compact, Cb_res, res_u, res_u_rowval, res_u_nzval)
            
            if Qa > Qmax
                Qmax = Qa
                aopt = compActs[aIdx]
            end
        end

      
        #Aopt[X_cidx] =  aopt yields a reference to compActs, which is BAD!
        nActsChanged = nActsChanged + ((Aopt[X_cidx][1] == aopt[1] && Aopt[X_cidx][2] == aopt[2]) ? 0 : 1)

        Aopt[X_cidx][1] = aopt[1]
        Aopt[X_cidx][2] = aopt[2]

        if (V_is_compact)
            maxVchange = max(maxVchange, abs(V[X_cidx] - Qmax))
            V[X_cidx] = Qmax
        else
            X_lidx = X2LIDX(X, Xp_indices[1])
            maxVchange = max(maxVchange, abs(V[X_lidx] - Qmax))
            V[X_lidx] = Qmax
            for i in 2:length(Xp_indices)
                X_lidx = X2LIDX(X, Xp_indices[i]) 
                V[X_lidx] = Qmax
            end
        end

        #TODO: remove this
#         if maxTime < Inf
#           elapsedTime = time() - start
#           if elapsedTime > maxTime
#             break
#           end
#         end
    end

    elapsedTime = time() - start
    if(maxVchange < 1 || iter==maxIters || elapsedTime > maxTime)
        @printf("Stopping after %i iterations (maxVchange = %.2f)\n", iter, maxVchange)
        break
    elseif mod(iter, 5) == 0
        @printf("At iteration #%i: maxVchange = %.2f, nActsChanged = %d, t = %.2f sec\n", 
                    iter, maxVchange, nActsChanged, elapsedTime)
    end
  end

  return Aopt

end


function greedyRandomRollout!(Qt_list, V::Vector{Float32}, ζ::Float32, β::Float32; maxIters::Int64=100, maxTime::Float64 = Inf)
    
    Aopt = (compActType)[copy(pattern.g_nullAct) for i in 1:g_nXcomp];
    
    V_is_compact = length(Aopt) == length(V)
    
    compActs = compActType[copy(pattern.g_nullAct) for i in 1:g_nMaxActs]
    
    Xp_indices = collect(permutations(1:g_nVehicles))
    
    n = size(Qt_list[1],1)
    res_u = spzeros(Float32, n,1)
    res_u_rowval = Array(Int64, n)
    res_u_nzval = Array(Float32, n)
    Cb_res = spzeros(Float32, n^(g_nVehicles-1),1); #n^K with K = nVehicles - 1
    
    
    start = time()
    @time for iter in 1:maxIters
      maxVchange = 0.0f0
      X_cidx = 0
      nActsChanged = 0;  
      for X in g_Xcomp
          X_cidx += 1 #X_cidx = X2CIDX(X)
          
          #Use the silent policy
          aopt = compActs[1]
          Va = QVeval(X, aopt, Qt_list, V, ζ, β, 
                      V_is_compact, Cb_res, res_u, res_u_rowval, res_u_nzval)
              
        
          if (V_is_compact)
              maxVchange = max(maxVchange, abs(V[X_cidx] - Va))
              V[X_cidx] = Va
          else
              X_lidx = X2LIDX(X, Xp_indices[1])
              maxVchange = max(maxVchange, abs(V[X_lidx] - Va))
              V[X_lidx] = Va
              for i in 2:length(Xp_indices)
                  X_lidx = X2LIDX(X, Xp_indices[i]) 
                  V[X_lidx] = Va
              end
          end
    
      end
    
      elapsedTime = time() - start
      if(maxVchange < 1 || iter==maxIters || elapsedTime > maxTime)
          @printf("Stopping after %i iterations (maxVchange = %.2f)\n", iter, maxVchange)
          break
      elseif mod(iter, 5) == 0
          @printf("At iteration #%i: maxVchange = %.2f, nActsChanged = %d, t = %.2f sec\n", 
                      iter, maxVchange, nActsChanged, elapsedTime)
      end
    end
    
    #Extract Aopt
    X_cidx = 0
    for X in g_Xcomp
         X_cidx += 1 #X_cidx = X2CIDX(X)
         
         aopt = pattern.g_nullAct
         Qmax = float32(-Inf)
         #Populate compact actions for this Xtate
         nActs = pattern.validCompActions!(compActs, X)
         for aIdx in 1:nActs
             Qa = QVeval(X, compActs[aIdx], Qt_list, V, ζ, β, 
                         V_is_compact, Cb_res, res_u, res_u_rowval, res_u_nzval)
             
             if Qa > Qmax
                 Qmax = Qa
                 aopt = compActs[aIdx]
             end
         end
    
         Aopt[X_cidx][1] = aopt[1]
         Aopt[X_cidx][2] = aopt[2]
     end

    return Aopt

end



function policy_X2a_compact(X::XType, Aopt::Vector{compActType})
  Xperm = sortperm(X)
  X_cidx = X2CIDX(X)

  compactAct = Aopt[X_cidx]

  act = pattern.g_nullAct
  if compactAct != pattern.g_nullAct
    pidx = [1:length(X)][Xperm[compactAct[1]]]
    act = [int8(pidx), compactAct[2]]
  end

  #This is still in compact form
  return act

end


function policy_X2a(X::XType, Aopt::Vector{compActType})
  #Get the compact form representation
  act = policy_X2a_compact(X, Aopt)
  #Trasnform it to the extended form
  return pattern.compAct2extAct(act,X2S(X))
end

function policy_S2a(S::SType, Aopt::Vector{compActType})
  return policy_X2a(S2X(S),Aopt)
end
##################################

Aopt = Array(compActType)
function savePolicy(Aopt, α, β_cost; prefix="")
    filename = "policies/" * prefix * "CTMDPpolicy_n_" * string(pattern.nPhases) * "_a_" * string(α) * "_b_" * string(β_cost) * ".jld"
    Aopt_idx = Int8[Aopt[i][1] for i in 1:length(Aopt)]
    Aopt_act = Int8[Aopt[i][2] for i in 1:length(Aopt)]
    JLD.save(filename, "Aopt_idx", Aopt_idx, "Aopt_act", Aopt_act); #"Vstar", Vshort,
    println(filename)
end

function loadPolicy(α, β_cost; prefix="")
    filename = "policies/" * prefix * "CTMDPpolicy_n_" * string(pattern.nPhases) * "_a_" * string(α) * "_b_" * string(β_cost) * ".jld"
    data = JLD.load(filename)
    global Aopt
    Aopt = compActType[[data["Aopt_idx"][i] , data["Aopt_act"][i]] for i in 1:length(data["Aopt_idx"])]
end

function ctmdpPolicy(S::Stype)
 return policy_S2a(S, Aopt::Vector{compActType})
end

##################################
#Use policy from the 1 phase case 
function liftUpPolicy!(Qt_list, V::Vector{Float32}, ζ::Float32, β::Float32,  Aopt1::Vector{compActType}; maxIters::Int64=100, maxTime::Float64 = Inf)
    #We are fine using references, as they will later get 
    #changed to other references, and since we ultimately populate
    #the value function, we won't return this to anyone!
    Aopt = Array(compActType, g_nXcomp);
    
    #We reconstruct the compact dictionary for phase free case!
    #TODO: Get rid of the hardcoded 27!!
    Scomp1 = combos_with_replacement(g_allstates[1:27], g_nVehicles)
    Xcomp1 = XType[S2X(s) for s in Scomp1]
    X2CIDX1_dict = (XType => Int64)[Xcomp1[cidx] => cidx for cidx in 1:length(Xcomp1)]
    function X1_2CIDX(X::XType)
      return X2CIDX1_dict[sort(X)]
    end

    X_cidx = 0
    for X in g_Xcomp
        X_cidx += 1 #X_cidx = X2CIDX(X)
        S = X2S(X)
        #Get rid of the phases since we are going to a 1phase case
        S1 = [phaseFree(s) for s in S]
        #Next convert to X1, note that this makes the assumption
        #that the phase free representation always comes first
        X1 = S2X(S1)
        #Next we need to figure out the compact index in the original
        #representation. Unfortunately we can't use X2CIDX since the dictionary
        #that is being used is the one with phases. So instead we use the locally
        #defined X1_2CIDX function
        X1_cidx = X1_2CIDX(X1)

        #Now we can go ahead and compute the optimal policy
        Aopt[X_cidx] = Aopt1[X1_cidx]
    end
    
    #Now use policy evaluation to figure out what V out-to-be
    Xp_indices = collect(permutations(1:g_nVehicles))
    n = size(Qt_list[1],1)
    res_u = spzeros(Float32, n,1)
    res_u_rowval = Array(Int64, n)
    res_u_nzval = Array(Float32, n)
    Cb_res = spzeros(Float32, n^(g_nVehicles-1),1); #n^K with K = nVehicles - 1
    
    
    start = time()
    @time for iter in 1:maxIters
        maxVchange = 0.0f0
        X_cidx = 0
        for X in g_Xcomp
            X_cidx += 1 #X_cidx = X2CIDX(X)
            Qa = QVeval(X, Aopt[X_cidx], Qt_list, V, ζ, β, false, 
                    Cb_res, res_u, res_u_rowval, res_u_nzval)
                
            X_lidx = X2LIDX(X, Xp_indices[1])
            maxVchange = max(maxVchange, abs(V[X_lidx] - Qa))
            V[X_lidx] = Qa
            for i in 2:length(Xp_indices)
                X_lidx = X2LIDX(X, Xp_indices[i]) 
                V[X_lidx] = Qa
            end
    
            #TODO: remove this
            if maxTime < Inf
              elapsedTime = time() - start
              if elapsedTime > maxTime
                break
              end
            end
        end
        
        elapsedTime = time() - start
        #Note that we use a much coarser value for maxVchange since
        #we just want to get close!
        if(maxVchange < 100 || iter==maxIters || elapsedTime > maxTime)
            @printf("Stopping after %i iterations (maxVchange = %.2f)\n", iter, maxVchange)
            break
        elseif mod(iter, 5) == 0
            @printf("At iteration #%i: maxVchange = %.2f , t = %.2f sec\n", 
                    iter, maxVchange, elapsedTime)
        end     
    end
end

###################

function uniformize(Q)
    T = copy(Q)
    dQ = diag(Q);
    κ = -minimum(dQ)
    Isp = spdiagm(ones(size(Q,1)))
    T = Q/κ + Isp
    #p = T * ones(size(T,1))
    #T = T * spdiagm(1./p) #normalize to 1
end
function sample(Q, dt)
    Isp = spdiagm(ones(size(Q,1)))
    M = Isp + Q * dt
end

end