
#It is assumed that someone else is providing the following:
#g_allStates: an array of symbols representing all possible substates that we can be in
#g_sn: a dictionary of type (Symbol => Int64) which goes from substate to index
#legalActions(S): function which tells us which actions are allowed!
#g_nullAct: the noaction case


#Number of nodes per instaces
g_nNodes = length(unique(g_allstates))

#Number of instances
g_nVehicles = 2


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

function X2LIDX(X::XType)
  #Note the reverse due to sub2ind using different convention than
  #what the kronecker sum results in!
  return sub2ind(g_XDIMS, reverse(X)...)
end
function S2LIDX(S::SType)
  return x2LIDX(S2X(S))
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
#Defining the Kronecker delta summation operator
function kronSum(A,B)
    b = size(B,1)
    a = size(A,2)
    Ib = spdiagm(ones(b))
    Ia = spdiagm(ones(a))

    return kron(A,Ib) + kron(Ia,B)
end
#Handling lists
function kronSum(Alist)
  A = Alist[1]
  for j = 2:length(Alist)
    A = kronSum(A, Alist[j])
  end
  return A
end


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


function QVeval(X::XType, action::typeof(g_nullAct), Qtjoint, Vcomp::Vector{Float64}, γ::Float64)

  β = 1./ γ

  (idx, act) = action;
  #Note the +1, we assume that act == 0 (from the null action) is a valid action!
  #i.e. all the Qtjoint arrays are offset by 1
  Qt = Qtjoint[act+1] #TODO: Make sure this is a reference and NOT a copy!


  #Next, re-order the states so that we are addressing the right instance
  Xo = swap(X, 1, idx)

  #Find the long index
  X_lidx  = X2LIDX(Xo)
  X_cidx = X2CIDX(Xo)



  #These are the rows of Qt for this column that are non-zero
  #which means these are the columns of Q for this
  #state (X).
  Qsparse_indices = Qt.colptr[X_lidx] : (Qt.colptr[X_lidx+1]-1)

  #Diagonal element is negative by construction
  #But q(x) is defined as the positive value
  qx = -Qt[X_lidx,X_lidx]

  #q(x) is not supposed to be part of the summation
  #we will later add Qt[si,si]*V[si] so starting with -Qt[si,si] will cancel it out
  qVsum =  qx * Vcomp[X_cidx]

  for Qsparse_idx in Qsparse_indices
    Qval = Qt.nzval[Qsparse_idx]
    V_LIDX = Qt.rowval[Qsparse_idx]
    V_CIDX = LIDX2CIDX(V_LIDX)

    qVsum += Qval * Vcomp[V_CIDX]
  end

  return (qVsum + r(X, action)) / (β + qx)
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
#############################################
function Reward(s::Vector{Symbol}, a::typeof(g_nullAct), β::Float64)
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

    (ncol, ntaxi) = NcolNtaxi(s)
    r += ncol * collisionCost + ntaxi * taxiCost;

    return r;
end

function r(X::XType, a::typeof(g_nullAct))
  return Reward(X2S(X), a, 0.01)
end

function gaussSeidel!(Vcomp::Vector{Float64}, γ::Float64)

  Aopt = (typeof(g_nullAct))[g_nullAct for i in 1:g_nXcomp];

  @time for iter in 1:100
    maxVchange = 0.
    for X in g_Xcomp
        aopt = g_nullAct
        Qmax = -Inf
        for a in legalActions(X2S(X))
            Qa = QVeval(X, a, Qt_joint, Vcomp, γ)
            if Qa > Qmax
                Qmax = Qa
                aopt = a
            end
        end

        X_cidx = X2CIDX(X)

        maxVchange = max(maxVchange, abs(Vcomp[X_cidx] - Qmax))
        Vcomp[X_cidx] = Qmax
        Aopt[X_cidx] = aopt
    end

    if(maxVchange < 1)
        @printf("Stopping after %i iterations (maxVchange = %.2f)\n", iter, maxVchange)
        break
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
