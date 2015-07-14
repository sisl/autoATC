module pattern

using Iterators

using auxFuns

export g_noaction, g_nullAct, isNullAct
export g_allstates, g_allstates_string 
export g_nNodes, g_nVehicles

export nPhases, phaseFree, phaseNum, phaseState, appendPhase
export randomChoice, NextStates

export s2x, S2X, x2s, X2S
export X2LIDX, S2LIDX, LIDX2X, LIDX2S

export xType, sType, XType, SType, compActType, extActType, ActType

rng = MersenneTwister()

#############################################
#Parameters
#############################################
const α = 1.0; #Probability of following ATC command

#Number of vehicle instances
const g_nVehicles = 4

#Number of phases per state
#TODO: Figure out why nPhases > 10 fails!
const nPhases = 1; #must be >= 1


###########################################
###########################################
const g_noaction = (int8(0), :∅)
const g_nullAct = Int8[0,0]

#Defining types
typealias sType Symbol; typealias SType Vector{sType}
typealias xType Int64 ; typealias XType Vector{xType}

typealias compActType typeof(pattern.g_nullAct)
typealias extActType  typeof(pattern.g_noaction)
typealias ActType Union(compActType,extActType)

function isNullAct(a::compActType)
    return a[1] == g_nullAct[1] #assumes anything with a[1] == 0 is null
end
function isNullAct(a::extActType)
    return a[1] == g_noaction[1] #assumes anything with a[1] == 0 is null
end



#############################################
#Set up states
const NextStates = (Symbol => Array{Symbol, 1})[]
#############################################

function addConn!(d, f, t)
    if haskey(d,f)
        push!(d[f], t)
        d[f] = unique(d[f])
    else
        d[f] = [t]
    end
end

normalTrans = { [:T,  :R, :U1, :LX1, :LD1, :LD2, :LB1, :F1, :R, :T];
                            [:U1, :RX1, :RD1, :RD2, :RB1, :F1];
                  [:F1, :GO, :U1, :U2, :LX2, :LD0, :LD1, :LD2, :LD3, :LB2, :F0, :F1];
                                 [:U2, :RX2, :RD0, :RD1, :RD2, :RD3, :RB2, :F0];
                                        [:U2, :LDep, :LArr, :LD1];[:LArr, :LD2];[:LArr, :LD3];
                                        [:U2, :RDep, :RArr, :RD1];[:RArr, :RD2];[:RArr, :RD3];}

allstates = [];
for k in 1:size(normalTrans, 1)
  for i in 2:length(normalTrans[k])
      addConn!(NextStates, normalTrans[k][i-1], normalTrans[k][i])
  end
  allstates = unique([allstates, normalTrans[k]])
end

#######################################################
#Add phases
#######################################################

#for now we will assume that we will use a fixed number of phases!
#Note that this will introduce (nPhases-1) as the last phase is
#assumed to be the state itself.
#We also assume that the actions are given at the last phase.
phaseFreeStates = [:R, :LDep, :LArr, :RDep, :RArr]
*(a::Symbol, b::Symbol) = symbol(string(a, b))
function appendPhase(s::Symbol, k::Int64)
    if isin(s, phaseFreeStates) || k >= nPhases || k <= 0
        return s
    else
        return symbol(string("ϕ",k,"_", s))
    end
end

function phaseState(s::Symbol)
    return string(s)[1] == 'ϕ'
end

const int0offset =  int('0')
function phaseNum(s::Symbol)
    if isin(s, phaseFreeStates)
        return 1
    else
        st = string(s)
        if st[1] == 'ϕ'
            return int(st[3])  - int0offset
        else
            return nPhases
        end
    end
end
function phaseFree(s::Symbol)
    st = string(s)
    if st[1] == 'ϕ'
        return symbol(st[5:end]) #This is going to bite you in the future!
    else
        return s
    end
end
if(nPhases > 1)
    for s in keys(NextStates)
        map!(s -> appendPhase(s,1), NextStates[s])
    end
    for s in allstates
        if !( s in phaseFreeStates)
            #Insert nPhases chain. Note that appendPhase(s,i+1)
            #returns s!
            for i in 1:(nPhases-1)
                NextStates[appendPhase(s,i)] = [appendPhase(s,i+1)]
            end
        end
    end
end
allstates = unique([allstates, collect(keys(NextStates))])
#######################################################
sn = (Symbol => Int64)[]
for i in 1:length(allstates)
    sn[allstates[i]] = i
end

const g_allstates = allstates;
const g_sn = sn;
const g_allstates_string = (UTF8String)[string(a) for a in g_allstates]

#Number of nodes per instaces
const g_nNodes = length(g_allstates)

###########################################
#Special states
###########################################
function listSpecialStates!(xOut, sOut, sList)
    cnt = 0
    for s in sList
      for k in 1:nPhases
        cnt += 1
        sOut[cnt] = appendPhase(s,k)
        xOut[cnt] = g_sn[sOut[cnt]]
      end
    end
    resize!(sOut, cnt)
    resize!(xOut, cnt)
end
const xDep = Array(xType, 2 * nPhases)
const sDep = Array(sType, 2 * nPhases)
listSpecialStates!(xDep, sDep, [:LDep, :RDep])
const xSafe = Array(xType, 5 * nPhases)
const sSafe = Array(sType, 5 * nPhases)
listSpecialStates!(xSafe, sSafe, [:LDep, :RDep, :LArr, :RArr, :T])
const xTaxi = Array(xType, 2 * nPhases)
const sTaxi = Array(sType, 2 * nPhases)
listSpecialStates!(xTaxi, sTaxi, [:T])
const xRunway = Array(xType, 2 * nPhases)
const sRunway = Array(sType, 2 * nPhases)
listSpecialStates!(xRunway, sRunway, [:R])
###########################################

###########################################
#Next states
###########################################
const NextXtates = Array(Vector{Int64}, length(g_allstates))
maxNextStates = 0
for x in 1:length(g_allstates)
  s = g_allstates[x]
  NextXtates[x] = [g_sn[sp] for sp in NextStates[s]]
  maxNextStates = max(maxNextStates, length(NextXtates[x]))
end
const g_nMaxActs = maxNextStates * g_nVehicles + 1
###########################################
#Add transition times for each state in minutes
###########################################

teaTime = (Symbol => Float64)[]

teaTime[:T]=230.26
teaTime[:R]=353.32
teaTime[:U1]=308.84
teaTime[:LX1]=138.13
teaTime[:LD1]=371.45
teaTime[:LD2]=332.55
teaTime[:LB1]=135.06
teaTime[:LX2]=136.00
teaTime[:LD0]=166.75
teaTime[:LD3]=227.07
teaTime[:LB2]=127.36
teaTime[:F0]=185.88
teaTime[:F1]=301.81
teaTime[:GO]=59.50
teaTime[:U2]=181.29
teaTime[:LDep]=588.25
teaTime[:LArr]=612.50

#Grab keys before we start inserting things
teaKeys = collect(keys(teaTime))
for s in teaKeys
    if !(s in phaseFreeStates)
        l = (teaTime[s] / nPhases);
        #Note that this will also modify teaTime[s]
        for i in 1:nPhases
            teaTime[appendPhase(s,i)] = l
        end
    end
end


# #Temporary hack, set all of them to the same value...
# for k in keys(teaTime)
#     teaTime[k] = 0.1;
# end

function symmetrize!(halfDict, symFun)
    for s in g_allstates
      if(string(phaseFree(s))[1] == 'R')
        astr = string(s)
        astr_l = replace(astr, 'R', 'L')
        a = symbol(astr)
        b = symbol(astr_l)
        if b in keys(halfDict)
          halfDict[a] = symFun(halfDict[b])
        end
      end
    end
end
symmetrize!(teaTime, x -> x)

#Make sure we didn't miss anything

assert(length(teaTime) == length(g_allstates))
#############################################
## Possible transitions
#############################################
function getNextPerms(ss::Array{Symbol,1})
  Nac = length(ss)
  ss_pos = [NextStates[s] for s in ss]

  ss_next = typeof(ss)[]
  for p in product(ss_pos...)
    push!(ss_next, Symbol[s for s in p])
  end

  return ss_next
end


#############################################
#coordinates
#############################################
dx = 2.; dy = 2.;
xy = Dict([:T, :R, :U1, :LX1, :LD1, :LD2, :LB1, :F1],
          {[0,-1] , [0, 0], [dx, 0], [1.5*dx, dy/2], [dx,dy], [-dx, dy], [-1.5*dx, dy/2], [-dx, 0]});
xy[:F0]  = [-2*dx, 0];
xy[:GO]  = [0, dy/2];
xy[:U2]  =  [2*dx, 0];
xy[:LX2] = [2.5*dx, dy/2];
xy[:LD0] = [2*dx, dy];
xy[:LD3] = [-2*dx, dy];
xy[:LB2] = [-2.5*dx, dy/2];
xy[:LDep] = [dx, 2*dy]
xy[:LArr] = [-dx, 2*dy]

symmetrize!(xy, x -> [x[1], -x[2]])

##
#############################################
##Printing to do file for vizualization
#############################################
#for f in keys(NextStates)
#    @printf("%s[label=\"%s(%i)\", pos=\"%.2f, %.2f\"]\n", f, f, sn[f], xy[f][1], xy[f][2])
#end
#for f in keys(NextStates)
#    for t in NextStates[f]
#        @printf("%s -> %s\n", f, t)
#    end
#end


#############################################
function randomChoice(from::Symbol, receivedATC::Bool, atcDesired::Symbol, rngState::AbstractRNG)
#############################################
#     snext = NextStates[from]
#     pnext = Float64[probFromTo(from, to, receivedATC, atcDesired) for to in snext]
# 
#     return snext[weightedChoice(pnext, rngState)]
    
    snext = NextStates[from]
    rnd = rand(rngState)
    
    N = length(snext)
    for k in 1:N
        rnd -= probFromTo(from, snext[k], receivedATC, atcDesired)
        if rnd < 0.0
            return snext[k]
        end
    end
    return snext[N]
    
end

#############################################
function probFromTo(from::Symbol, to::Symbol, receivedATC::Bool, atcDesired::Symbol)
#############################################
    p = 0.
    allNext = NextStates[from]

    if isin(to, allNext) 
        Nnext = length(allNext);
        #Only one option, use it
        if Nnext == 1
            p = 1.
        #Not addressing this aircraft or invalid atc command
        #All other states are equally likely
        elseif !receivedATC || !(atcDesired in allNext)
            p = 1./ Nnext
                        
            #Ugly handling of Taxi/Runway interaction!           
            if isin(from, pattern.sRunway)
                if isin(to, pattern.sTaxi)
                    p = 0.1
                else
                    p = 0.9
                end
            elseif isin(from, pattern.sTaxi)
                if isin(to, pattern.sTaxi)
                    p = 0.1
                else
                    p = 0.9
                end
            end
        #received ATC command, collaborate!
        else
          if to == atcDesired
                p = α::Float64
          else
                p = (1-α::Float64)/(Nnext-1)
          end
        end
    end
    return p
end




#############################################
function validActions(S)
#############################################
  A = [g_noaction]; sizehint(A, 10);
  for (i, s) in enumerate(S)
    #Can't tell departing aircrafts what to do
    if !(s in [:LDep, :RDep, :R])
      snext = NextStates[s]
      if(length(snext) > 1)
        for sn in snext
          #Can't tell aircraft to depart...
          if !(sn in [:LDep, :RDep])
            push!(A, (i, sn))
          end
        end
      end
    end
  end
  return A
end


function validActions!(acts::Vector{typeof(g_noaction)}, S::Vector{Symbol})
  nActs = 0

  nActs += 1
  #acts[nActs] = copy(g_noaction) #should be the case by default
  for i in 1:length(S)
    s = S[i]
    #on runway or in departure state, can't tell it what to do
    if isin(s, sDep) || isin(s, sRunway)
      continue
    end
    
    snext = NextStates[s]
    nX = length(NextStates[s])
    
    if nX > 1
      for j in 1:nX
        if isin(snext[j], sDep) #Can't tell aircraft to depart
          continue
        end
        nActs += 1
        acts[nActs] = (i, snext[j])
      end
    end
  end
  return nActs
end


function extAct2compAct(act::typeof(g_noaction), S)
  cAct = copy(g_nullAct)
  if act != g_noaction
    sidx = act[1]
    cAct[1] = sidx
    for (cidx, s) in enumerate(NextStates[S[sidx]])
      if s == act[2]
        cAct[2] = cidx
      end
    end
  end
  return cAct
end

function compAct2extAct(act::typeof(g_nullAct), S)
  eAct = copy(g_noaction)
  if act != g_nullAct && act[1] > 0 && act[1] <= length(S)
    sidx = act[1]
    Sp = NextStates[S[sidx]]
    if act[2] <= length(Sp) #otherwise not a valid action!
      eAct = (act[1], Sp[act[2]])
    end
  end
  return eAct
end

function validCompactActions(S)
  actions = validActions(S)
  compActions = Array(typeof(g_nullAct),length(actions),1)
  for (idx, act) in enumerate(actions)
    compActions[idx] = extAct2compAct(act, S)
  end
  return compActions
end

function validCompActions!(compActs::Vector{typeof(g_nullAct)}, X::Vector{Int64})
  nActs = 0

  nActs += 1
  #compActs[nActs] = copy(g_nullAct) #should be the case by default
  for i in 1:length(X)
    x = X[i]
    #on runway or in departure state, can't tell it what to do
    #Note that if a similar x has been encountered, no point in evaluating
    #commanding the same one
    #x in X[1:(i-1)] #This causes memory being allocated ...
    #Instead using little for loop
    isrepeated = false
    for j in 1:(i-1)
        (isrepeated = (x == X[j])) && break
    end
    if isrepeated || isin(x, xDep) || isin(x, xRunway) #||  x in X[1:(i-1)]
      continue
    end
    nX = length(NextXtates[x])
    if nX > 1
      for j in 1:nX
        if isin(NextXtates[x][j], xDep) #Can't tell aircraft to depart
          continue
        end
        nActs += 1
        compActs[nActs][1] = i
        compActs[nActs][2] = j
      end
    end
  end
  return nActs
end
#############################################
function Transition(S::Array{Symbol,1}, a::typeof(g_noaction), Snext::Array{Symbol,1})
#############################################
    p = 1.;
    idx = 1;
    for (from, to) in zip(S, Snext)
        p *= probFromTo(from, to, idx == a[1], a[2])
        idx += 1;
    end
    return p;
end

function Transition(S::Array{Symbol,1}, acomp::typeof(g_nullAct), Snext::Array{Symbol,1})
  aext = compAct2extAct(acomp, S)
  return Transition(S, aext, Snext)
end

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

#Going between symbol representation and index representation
function s2x(s::sType)
  return (g_sn::Dict{sType, xType})[s]
end
function S2X(S::SType)
  return xType[(g_sn::Dict{sType,xType})[s] for s in S]
end
function x2s(x::xType)
  return (g_allstates::Vector{sType})[x]
end
function X2S(X::XType)
  return sType[(g_allstates::Vector{sType})[x] for x in X]
end

#Long indices are obtained from sub2ind/ind2sub
#where the dimensions are governed by the number
#of nodes and number of vehicles!
const g_XDIMS = tuple((g_nNodes*ones(Int64, g_nVehicles))...)

function X2LIDX(X::XType)
    #Note the reverse due to sub2ind using column major numbering whereas
    #the kronecker math uses a row major numbering
    index = X[g_nVehicles]
    stride = 1
    for k in (g_nVehicles-1):-1:1
        stride = stride * g_nNodes
        index += (X[k]-1) * stride
    end
    return index
end

#same as above, but takes permutation into account
function X2LIDX(X::XType, perm::XType)
    index = X[perm[g_nVehicles]]
    stride = 1
    for k = (g_nVehicles-1):-1:1
        stride = stride * g_nNodes
        index += (X[perm[k]] -1) * stride
    end
    return index
end

function S2LIDX(S::SType)
  return X2LIDX(S2X(S))
end

function LIDX2X(lidx::Int64)
  #TODO: make this faster ...
  X_rev = xType[ind2sub(g_XDIMS, lidx) ...]
  return reverse(X_rev)
end
function LIDX2S(lidx::Int64)
  X = LIDX2X(lidx)
  return X2S(X)
end

#############################################
#Old simulation code from MDP...
#TODO: Delete this?
#############################################
function simulate(s; policy::Dict = Dict(), N=10)
  return simulate(s, (s-> policy[s]), N)
end

function simulate(s, policy::Function; N=10, endEarly = false)
#############################################
    S = {s}
    A = {g_noaction}
    for step in 1:N
        Snow = S[end]
        a = policy(Snow)
        Snew = deepcopy(Snow)
        for i in 1:length(Snow)
            Snew[i] = randomChoice(Snow[i], a[1] == i, a[2], rng)
        end
        push!(S, Snew)
        push!(A, a)
        if(endEarly && NcolNtaxi(Snew)[1] != 0)
          break;
        end
    end
    A = A[2:end];
    rewards = Float64[Reward(s, a) for (s, a) in zip(S,A)];
    collisions = [NcolNtaxi(s)[1] for s in S]
    return S, A, collisions, rewards
end



#############################################
function plotSim(S, Actions, collisions, rewards; animation = false, Ncolprev = 0)
#############################################
    M = "";

    xrange = 6.; yrange = 5.;
    if(animation)
      M = "o"
      ax = PyPlot.subplot(1,1,1)
    else
      ax = PyPlot.subplot(2,1,1)
    end

    for s in g_allstates
      (x, y) = xy[s]
        PyPlot.text(x, y, s, fontsize=10, alpha=0.8, bbox={"facecolor" => "white", "alpha"=>0.2})
    end

    colors = ["blue", "green", "black", "red"]
    for n in 1:length(S[1])
        x = [xy[s[n]][1] for s in S]+(n-1)*.1
        y = [xy[s[n]][2] for s in S]+(n-1)*.05
        PyPlot.plot(x, y, linestyle="--", linewidth = 2, color=colors[n]);
        if(animation)
          M = (n == Actions[end][1]) ? "D" : "o"
          PyPlot.plot(x[end], y[end], linestyle="", color=colors[n], marker=M);
          if(M == "D")
             text(x[end], y[end], Actions[end][2], fontsize=12)
          end
        end
    end


    cidx = find(x -> x >= 1, collisions)
    Nc = float(length(cidx));
    NcolSum = cumsum(collisions);
    RewardSum = cumsum(rewards)/1000.;

    NcCount = Dict(g_allstates, zeros(size(g_allstates)));
    if(Nc > 0.)
      if(animation)
        cidx = [cidx[end]]
      end
      Nsum = 0
      for ss in S[cidx]
        for s in unique(ss)
          if s != :T
            NcCount[s] += (c = count(x -> x == s, ss) - 1)
            Nsum += c
          end
        end
      end

      for s in g_allstates
            msize = int(50.*NcCount[s]/Nsum)
            x = xy[s][1]; y = xy[s][2];
            if animation
                msize = min(msize, 10)
            end
            if msize > 0
                PyPlot.plot(x, y, marker="o", color="red", markersize = msize, alpha = 0.9);
            end
      end
    end

    if animation
      #PyPlot.title("Total Collisions = " * string(collisions[end] + Ncolprev))
    end
    ax[:set_xlim]([-.9, 1]*xrange);
    ax[:set_ylim]([-1, 1]*yrange);
    ax[:set_yticklabels]([]); ax[:set_xticklabels]([])
    #PyPlot.grid("on")



    if (!animation)
      ax = PyPlot.subplot(2,1,2)
      PyPlot.plot(-NcolSum, color="red")
      PyPlot.plot(RewardSum, color="blue")
      ax[:set_ylim]([minimum([RewardSum, NcolSum]), 0.01]);
      PyPlot.grid("on")
      PyPlot.legend(["Number of Collisions", "R / 1000."])
    end

end



end

