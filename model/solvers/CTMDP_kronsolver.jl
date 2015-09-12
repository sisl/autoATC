module CTMDP_kronsolver

using pattern
using RewardFun
using auxFuns
using kronfun

using CTMDP_kronindexing
using CTMDP_kron

export solveCTMDP

function QVeval(X::XType, action::compActType, Qt_list, V::Vector{Float32},
                ζ::Float32, β::Float32, V_is_compact::Bool, 
                Cb_res, res_u, res_u_rowval, res_u_nzval)

  R = Reward(X, action, β)

  #This is a terminal state... Q(s,a) = R(s,a)
  assert(β < 0.9f0) #We make the assumption that action cost is small relative to collision cost
  if( R <=  RewardFun.collisionCost)
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
  Qt = Qti_ABt!(Qt_list[act+1], Qt_list[1], g_nVehicles-1, X_lidx, 
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
    nActsChanged = 0;
    X = zeros(xType, g_nVehicles)
    for X_cidx in 1:g_nXcomp
        CIDX2X!(X, X_cidx)
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
      X = zeros(xType, g_nVehicles)
      for X_cidx in 1:g_nXcomp
          CIDX2X!(X, X_cidx)
              
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


####################
#Putting it all together...
#[0.0f0, 0.001f0, 0.005f0, 0.01f0, 0.05f0]
function solveCTMDP(β_costs=[0.0f0, 0.001f0, 0.005f0, 0.01f0, 0.05f0],
                    ζ_discount=float32(0.5/60.))

    amax = maximum([length(pattern.NextStates[s]) for s in keys(pattern.NextStates)]);
    A = (Int8)[0:amax];
    P0 = spzeros(Float32, pattern.g_nNodes, pattern.g_nNodes);
    P = (typeof(P0))[copy(P0) for a in A]
    for a in A
        if a == 0
            act = pattern.g_nullAct
        else
            act = Int8[1, a]
        end
        for x in 1:pattern.g_nNodes
            s = x2s(x)
            for sp in pattern.NextStates[s]
                xp = s2x(sp)
                P[a+1][x,xp] = pattern.Transition([s], act, [sp])
            end
        end
    end
    M0 = speye(Float32, g_nNodes)
    for x in 1:g_nNodes
        s = x2s(x)
        M0[x,x] = 1./(pattern.sojurnTime[s])
    end
    
    Isp = speye(Float32, g_nNodes);
    Qt_list = (typeof(M0))[(M0*(P[a+1] - Isp))' for a in A];
    Vlong = zeros(Float32,g_nSlong);

    for β_cost in β_costs 
        println("Solving for ", β_cost)
        Aopt = gaussSeidel!(Qt_list, Vlong, ζ_discount, β_cost); #maxIters = 1, maxTime = 5.*30.);
        saveCTMDPpolicy(Aopt, pattern.α, β_cost)
    end

end



end