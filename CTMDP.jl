include("kronfun.jl")


#It is assumed that someone else is providing the following:
#g_allStates: an array of symbols representing all possible substates that we can be in
#g_sn: a dictionary of type (Symbol => Int64) which goes from substate to index
#legalActions(S): function which tells us which actions are allowed!
#g_nullAct: the noaction case


#Number of nodes per instaces
const g_nNodes = length(unique(g_allstates))

#Number of instances
const g_nVehicles = 4


function combos_with_replacement(list, k)
    n = length(list)
    [[list[c[i]-i+1] for i=1:length(c)] for c in combinations([1:(n+k-1)],k)]
end


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

#Defining types
const sType = Symbol; const SType = Vector{Symbol}
const xType = Int64; const XType = Vector{Int64}

#Going between symbol representation and index representation
function s2x(s::sType)
  return (g_sn::Dict{Symbol, Int64})[s]
end
function S2X(S::SType)
  return xType[(g_sn::Dict{Symbol,Int64})[s] for s in S]
end
function x2s(x::xType)
  return (g_allstates::Vector{Symbol})[x]
end
function X2S(X::XType)
  return sType[(g_allstates::Vector{Symbol})[x] for x in X]
end


#Get all possible states, allowing replacement, but order does not matter:
const g_Scomp = combos_with_replacement(g_allstates, g_nVehicles)
const g_Xcomp = XType[S2X(s) for s in g_Scomp]

const g_nScomp = length(g_Scomp); const g_nXcomp = g_nScomp;
const g_nSlong  = g_nNodes^g_nVehicles; const g_nXlong = g_nSlong;

const g_nCompActs = maxNextStates * g_nVehicles + 1

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

#Long indices are obtained from sub2ind/ind2sub
#where the dimensions are governed by the number
#of nodes and number of vehicles!


const g_XDIMS = tuple((g_nNodes*ones(Int64, g_nVehicles))...)


function Xsub2ind(X::XType)
    index = X[1]
    stride = 1
    for k=2:g_nVehicles
        stride = stride * g_nNodes
        index += (X[k]-1) * stride
    end
    return index
end

function X2LIDX(X::XType)
  #Note the reverse due to sub2ind using different convention than
  #what the kronecker sum results in!
  return Xsub2ind(reverse(X))
end

#same as above, but takes permutation into account
function X2LIDX(X::XType, perm::XType)
    index = X[perm[g_nVehicles]]
    stride = 1
    for k= (g_nVehicles-1):-1:1
        stride = stride * g_nNodes
        index += (X[perm[k]] -1) * stride
    end
    return index
end


function S2LIDX(S::SType)
  return X2LIDX(S2X(S))
end

function LIDX2X(lidx::Int64)
  X_rev = xType[ind2sub(g_XDIMS, lidx) ...]
  return reverse(X_rev)
end
function LIDX2S(lidx::Int64)
  X = LIDX2X(lidx)
  return X2S(X)
end


#This is when we need to go between compact and long indices
#Note that many LIDX map to the same CIDX
#And that one CIDX maps to many LIDX (but we only return one)
function LIDX2CIDX(lidx::Int64)
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


###########################################
#Generic functions
###########################################


function swap{T}(s::Vector{T}, i, j)
  sc = copy(s);
  if(i != 0 && j != 0)
    sc[i] = s[j];
    sc[j] = s[i];
  end

  return sc
end

###########################################


function findn_rows{Tv,Ti}(S::SparseMatrixCSC{Tv,Ti}, colIdx::Integer)
    idx = S.colptr[colIdx] : (S.colptr[colIdx+1]-1)
    return (S.rowval[idx] , S.nzval[idx])
end


function QVeval(X::XType, action::typeof(g_nullAct), Qt_list, V::Vector{Float32}, β::Float64, V_is_compact::Bool, res_u_rowval, res_u_nzval)

  (idx, act) = action;

  #First, re-order the states so that we are addressing the right instance
  Xo = swap(X, 1, idx)

  #Find the long index
  X_cidx = X_lidx  = X2LIDX(Xo)
  if(V_is_compact)
    X_cidx = X2CIDX(Xo)
  end 
  
  
  #get the relevant Q row
  Qt = Qti_ABt(Qt_list[act+1], Qt_list[1], g_nVehicles-1, X_lidx, res_u_rowval, res_u_nzval)

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

  return (qVsum ) / (β + qx) + r(X, action)

end


#############################################
function NcolNtaxi(s::Vector{Symbol})
#############################################
  #We don't need to worry about collisions
  #if we are in taxi, departure or arrival
  notColl = find(x -> !(x in [:T, :LDep, :RDep, :LArr, :RArr]), s)
  inTaxi = find(x -> x == :T, s)

  Nc = length(notColl);
  Ncu = length(unique(s[notColl]));

  return [Nc - Ncu, length(inTaxi)]
end

function NcolNtaxi(X::XType)
  Nc = 0
  Nt = 0
  for i in 1:length(X)
    if X[i] in xTaxi
      Nt += 1
    end
    if X[i] in xSafe
      continue
    end
    for j in (i+1):length(X)
      if X[j] == X[i]
        Nc += 1
      end
    end
  end
  return [Nc, Nt]
end
#############################################
function Reward(s::Vector{Symbol}, a::typeof(g_nullAct), β::Float64)
  return Reward(S2X(S), a, β)
end
function Reward(X::XType, a::typeof(g_nullAct), β::Float64)
#############################################
    r = 0.

    #Each collision costs 1000.
    #And each aircraft just sitting on the taxi also incurs cost
    collisionCost = -1000.
    taxiCost = -10.

    #Actions have a cost
    if(a != g_nullAct)
        r += β * collisionCost;
    end

    (ncol, ntaxi) = NcolNtaxi(X)
    r += ncol * collisionCost + ntaxi * taxiCost;

    return r;
end

function r(X::XType, a::typeof(g_nullAct))
  return Reward(X, a, β_cost)
end

function gaussSeidel!(Qt_list, V::Vector{Float32}, β::Float64; maxIters::Int64=100, maxTime::Float64 = Inf)

  Aopt = (typeof(g_nullAct))[g_nullAct for i in 1:g_nXcomp];
  
  V_is_compact = length(Aopt) == length(V)

  compActs = Array(typeof(g_nullAct), g_nCompActs)

  Xp_indices = collect(permutations(1:g_nVehicles))
  
  n = size(Qt_list[1],1)
  res_u_rowval = Array(Int64, n)
  res_u_nzval = Array(Float64, n)

  start = time()
  @time for iter in 1:maxIters
    maxVchange = 0.
    for X in g_Xcomp
        aopt = g_nullAct
        Qmax = -Inf
        #Populate compact actions for this Xtate
        nActs = validCompActions!(compActs, X)
        for aIdx in 1:nActs  #in legalActions(X2S(X))
            Qa = QVeval(X, compActs[aIdx], Qt_list, V, β, V_is_compact, res_u_rowval, res_u_nzval)
            if Qa > Qmax
                Qmax = Qa
                aopt = compActs[aIdx]
            end
        end

        X_cidx = X2CIDX(X)
        Aopt[X_cidx] = aopt
    
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
        if maxTime < Inf
          elapsedTime = time() - start
          if elapsedTime > maxTime
            break
          end
        end
    end

    elapsedTime = time() - start
    if(maxVchange < 1 || iter==maxIters || elapsedTime > maxTime)
        @printf("Stopping after %i iterations (maxVchange = %.2f)\n", iter, maxVchange)
        break
    elseif mod(iter, 5) == 0
        @printf("At iteration #%i (%.2f sec)\n", iter, elapsedTime)
    end
  end

  return Aopt

end


function policy_X2a(X::XType, Aopt::Vector{typeof(g_nullAct)})
  Xperm = sortperm(X)
  X_cidx = X2CIDX(X)

  compactAct = Aopt[X_cidx]

  act = g_nullAct
  if compactAct != g_nullAct
    pidx = [1:length(X)][Xperm[compactAct[1]]]
    act = [pidx, compactAct[2]]
  end

  try
    act = compAct2extAct(act,X2S(X))
  catch
    println(X, sort(X), act, compactAct)
  end
  return act

end
function policy_S2a(S::SType, Aopt::Vector{typeof(g_nullAct)})
  return policy_X2a(S2X(S),Aopt)
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
