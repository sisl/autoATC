using HDF5, JLD

if isdefined(current_module(),:LastMain) && isdefined(LastMain,:Iterators)
  using LastMain.Iterators
else
  using Iterators
end




#############################################
#Parameters
#############################################
α = 0.95; #Probability of following ATC command
β = 0.01; #Fraction of atc cost relative to collision cost
noaction = (0, :∅)

#############################################
#Set up states
NextStates = (Symbol => Array{Symbol, 1})[]
#############################################

function addConn!(d, f, t)
    if haskey(d,f)
        push!(d[f], t)
        d[f] = unique(d[f])
    else
        d[f] = [t]
    end
end

normalTrans = { [:T, :T, :R, :U1, :LX1, :LD1, :LD2, :LB1, :F1, :R, :T];
                            [:U1, :RX1, :RD1, :RD2, :RB1, :F1];
                  [:F1, :GO, :U1, :U2, :LX2, :LD0, :LD1, :LD2, :LD3, :LB2, :F0, :F1];
                                 [:U2, :RX2, :RD0, :RD1, :RD2, :RD3, :RB2, :F0];
                                        [:U2, :LDep, :LDep, :LArr, :LArr, :LD1];[:LArr, :LD2];[:LArr, :LD3];
                                        [:U2, :RDep, :RDep, :RArr, :RArr, :RD1];[:RArr, :RD2];[:RArr, :RD3];}

allstates = [];
for k in 1:size(normalTrans, 1)
  for i in 2:length(normalTrans[k])
      addConn!(NextStates, normalTrans[k][i-1], normalTrans[k][i])
  end
  allstates = unique([allstates, normalTrans[k]])
end
sn = (Symbol => Int64)[]
for i in 1:length(allstates)
    sn[allstates[i]] = i
end



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
allstates_string = [string(a) for a in allstates]
for astr in allstates_string
  if(astr[1] == 'R')
    astr_l = replace(astr, 'R', 'L')
    a = symbol(astr)
    b = symbol(astr_l)
    if b in keys(xy)
      xy[a] = [xy[b][1], -xy[b][2]]
    end
  end
end


#############################################
#Probability of following ATC command
#############################################
function weightedChoice(weights)
    rnd = rand() * sum(weights)
    for (i, w) in enumerate(weights)
        rnd -= w
        if rnd < 0
            return i
        end
    end
end

#############################################
function randomChoice(from::Symbol, receivedATC::Bool, atcDesired::Symbol)
#############################################
    snext = NextStates[from]
    pnext = Float64[probFromTo(from, to, receivedATC, atcDesired) for to in snext]

    return snext[weightedChoice(pnext)]
end

#############################################
function probFromTo(from::Symbol, to::Symbol, receivedATC::Bool, atcDesired::Symbol)
#############################################
    p = 0.
    allNext = NextStates[from]
    if to in allNext
        Nnext = length(allNext);
        #Only one option, use it
        if Nnext == 1
            p = 1.
        #Not addressing this aircraft or invalid atc command
        elseif !receivedATC || !(atcDesired in allNext)
            #except for taxi (aircraft will not take-off unless told so)
            #All other states are equally likely
            if(from != :T)
              p = 1./ Nnext
            #i.e. from Taxi -> Taxi without ATC commands
            elseif (to == :T)
                p = α
            #From Taxi -> other things..
            else
                p = (1-α)/(Nnext-1)
            end
        #received ATC command, collaborate!
        else
          if to == atcDesired
            p = α
          else
            p = (1-α)/(Nnext-1)
          end
        end
    end
    return p
end


#############################################
function simulate(s; policy::Dict = Dict(), N=10)
  return simulate(s, (s-> policy[s]), N)
end

function simulate(s, policy::Function; N=10, endEarly = false)
#############################################
    S = {s}
    A = {noaction}
    for step in 1:N
        Snow = S[end]
        a = policy(Snow)
        Snew = deepcopy(Snow)
        for i in 1:length(Snow)
            Snew[i] = randomChoice(Snow[i], a[1] == i, a[2])
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

    for s in allstates
      (x, y) = xy[s]
      PyPlot.text(x, y, s, fontsize=8, alpha=0.5)
    end

    colors = ["blue", "green", "black", "orange"]
    for n in 1:length(S[1])
        x = [xy[s[n]][1] for s in S]+(n-1)*.1
        y = [xy[s[n]][2] for s in S]
        PyPlot.plot(x, y, linestyle="--", color=colors[n]);
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

    NcCount = Dict(allstates, zeros(size(allstates)));
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

      for s in allstates
            msize = int(50.*NcCount[s]/Nsum)
            x = xy[s][1]; y = xy[s][2];
            if animation
                msize = min(msize, 10)
            end
            if msize > 0
                PyPlot.plot(x, y, marker="o", color="red", markersize = msize, alpha = 0.8);
            end
      end
    end

    if animation
      PyPlot.title("Total Collisions = " * string(collisions[end] + Ncolprev))
    end
    ax[:set_xlim]([-1, 1]*xrange);
    ax[:set_ylim]([-1, 1]*yrange);
    PyPlot.grid("on")



    if (!animation)
      ax = PyPlot.subplot(2,1,2)
      PyPlot.plot(-NcolSum, color="red")
      PyPlot.plot(RewardSum, color="blue")
      ax[:set_ylim]([minimum([RewardSum, NcolSum]), 0.01]);
      PyPlot.grid("on")
      PyPlot.legend(["Number of Collisions", "R / 1000."])
    end

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
function Reward(s::Vector{Symbol}, a::typeof(noaction))
#############################################
    r = 0.

    #Each collision costs 1000.
    #And each aircraft just sitting on the taxi also incurs cost
    collisionCost = -1000.
    taxiCost = -10.

    #Actions have a cost
    if(a != noaction)
        r += β * collisionCost;
    end

    (ncol, ntaxi) = NcolNtaxi(s)
    r += ncol * collisionCost + ntaxi * taxiCost;

    return r;
end



#############################################
function validActions(ss)
#############################################
  A = [noaction]; sizehint(A, 10);
  for (i, s) in enumerate(ss)
    #Can't tell departing aircrafts what to do
    if !(s in [:LDep, :RDep])
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


#############################################
function Transition(s::Array{Symbol,1}, a::typeof(noaction), snext::Array{Symbol,1})
#############################################
    p = 1.;
    idx = 1;
    for (from, to) in zip(s, snext)
        p *= probFromTo(from, to, idx == a[1], a[2])
        idx += 1;
    end
    return p;
end

#############################################
#Defining states
#############################################


if(false)

@printf("Defining states \n"); tic()
S = Vector{Symbol}[[a,b,c,d] for a in allstates, b in allstates, c in allstates, d in allstates]; S = S[:];
NS = length(S)
s2i = Dict(S,[1:NS])

toc(); @printf("Defining permutations/actions \n"); tic();
NextPerms = typeof(getNextPerms(S[1]))[getNextPerms(s) for s in S]
validActs = typeof(validActions(S[1]))[validActions(s) for s in S]
numValidActs = Int8[length(validActs[i]) for i in 1:NS]


toc(); @printf("Allocating memory for transitions/Rewards \n"); tic();
tsaEntry = typeof(Float32[])[]
rsaEntry = Float32[];
TSA = Array(typeof(tsaEntry), NS)
RSA = Array(typeof(rsaEntry), NS)
Jlist = Array(typeof(Int64[]), NS)
@time for from_i in 1:NS
  if(from_i % int(NS/4) == 0)
    toc(); @printf("%i/%i ... \n",from_i, NS); tic();
  end

  from = S[from_i]
  nn =  length(validActs[from_i])

  RSA[from_i] = Array(Float32, nn)
  TSA[from_i] = Array(typeof(Float32[]), nn)

  from_i_nextPerms = NextPerms[from_i]
  Jlist[from_i] = Int64[s2i[to] for to in from_i_nextPerms]

  for (a_i, a) in enumerate(validActs[from_i])
      TSA[from_i][a_i] = Array(Float32, length(from_i_nextPerms))
  end

end
toc();

############################################
#Get a better ordering
#############################################
# function breadthFirst!(s0::Vector{Symbol}, Ordering::Vector{Int32})
#     breadthFirst!(int32(s2i[s0]), Ordering)
# end
# function breadthFirst!(s0_i::Int32, Ordering::Vector{Int32})
#     visited = falses(NS)
#     stack = Int32[]
#     function bless(i::Int32)
#         push!(stack, i)
#         push!(Ordering, i)
#         visited[i] = true
#     end
#     bless(s0_i)
#     while(length(stack) > 0 && length(Ordering) < NS)
#         scurrent_i = stack[end]
#         nochildren = true
#         for s1 in NextPerms[scurrent_i]
#             s1_i = int32(s2i[s1])
#             if !visited[s1_i]
#                 bless(s1_i)
#             end
#         end
#         pop!(stack)
#     end
# end

# @printf("defining breadth first ordering \n")
# o = Int32[]; sizehint(o, NS)
# @time breadthFirst!(S[1], o)
# reverse!(o)
o = Int32[1:NS]
#############################################
#Updating TSA and RSA
#############################################
function UpdateTSA!(TSA::Vector{typeof(tsaEntry)}, validActs::Vector{Vector{typeof(noaction)}}, NextPerms::Vector{Vector{Vector{Symbol}}})
  for from_i in 1:NS
    if(from_i % int(NS/5) == 0)
      @printf("%i/%i ... ",from_i, NS);
    end

    from_i_nextPerms = NextPerms[from_i]
    for a_i in 1:length(validActs[from_i])
      a = validActs[from_i][a_i]
      for to_i in 1:length(NextPerms[from_i])
        to = from_i_nextPerms[to_i]
        TSA[from_i][a_i][to_i] = Transition(S[from_i], a, to)
      end
    end
  end
end


function UpdateRSA!(RSA::Vector{typeof(rsaEntry)}, validActs::Vector{Vector{typeof(noaction)}})
  for from_i in 1:NS
    if(from_i % int(NS/5) == 0)
      @printf("%i/%i ... ",from_i, NS);
    end

    from = S[from_i]::Vector{Symbol}

    for a_i in 1:length(validActs[from_i])
      RSA[from_i][a_i] = Reward(from, validActs[from_i][a_i])
    end
  end
end

#############################################
#Dot product for sparse matrices (sort of)
#############################################
function sparsedot(U::Array{Float32,1}, M::Array{Float32,1}, s0_i::Int32)
    s = 0.0f0
    i = 1
    for idx in Jlist[s0_i]::Array{Int64, 1} #(i, m) in zip(M.rowval, M.nzval)
        @inbounds s += U[idx] * M[i]
        i += 1
    end
    return s
end


#############################################
#Gauss Seidel Value Iteration
#############################################
function gaussSeidelValueIteration(numIterations::Integer, Uarray::Array{Float32,1};  γ = 0.95f0)
    Piarray = Array(Int64, NS)

    function Q(s0_i::Int32, a_i::Int8)
        @inbounds return RSA[s0_i][a_i]::Float32 + γ * sparsedot(Uarray, TSA[s0_i][a_i], s0_i::Int32)::Float32
    end

    for t in 1:numIterations
        maxDelta = 0.0f0;
        for s0_i::Int32 in o
            amax = 1; Qmax = -Inf32;
            for a_i::Int8 in 1:numValidActs[s0_i]::Int8
                Qnew = Q(s0_i, a_i)::Float32
                if Qnew > Qmax
                    Qmax = Qnew
                    amax = a_i
                end
            end
            @inbounds maxDelta = max(maxDelta, abs(Qmax - Uarray[s0_i]))
            @inbounds Uarray[s0_i] = Qmax;
            @inbounds Piarray[s0_i] = amax;
        end

        if  maxDelta < 1
            @printf("No significant change in U after %i iterations. Terminating\n", t)
            break
        elseif t == numIterations
            @printf("Ran all %i iterations and maxDelta == %.4f\n", t, maxDelta)
        elseif (t % 5) == 1
            @printf("iter #%i ... maxDelta = %.4f\n", t, maxDelta)
        end
    end
    return Uarray, Piarray

end


#############################################
#Policy function
#############################################
function PiStarFun(s0::Array{Symbol, 1})
    s_i = s2i[s0];
    return validActs[s_i][Pistar[s_i]]
end

#############################################
#Solving
#############################################
if !isdefined(current_module(),:solveAll)
  solveAll = false
end

α_range = [0., 0.05, 0.25, 0.5, 0.75, 0.85, 0.95, 1.]
β_range = [0., 0.01, 0.1, 0.25, 0.5, 1.]

Ustar = zeros(Float32, NS);

tic();
if(solveAll)
  for i in 1:length(α_range)
    α = α_range[i];
    @printf("Updating TSA for α=%.2f: ", α);
    @time UpdateTSA!(TSA, validActs, NextPerms);

    for j in 1:length(β_range)
      β = β_range[j];
      @printf("Updating RSA for β=%.2f: ", β);
      @time UpdateRSA!(RSA, validActs);

      @printf("α = %.2f, β = %.2f\n", α, β)
      #Hot starting ...
      @printf("Running Gauss: \n");
      @time (Ustar, Pistar) = gaussSeidelValueIteration(50, Ustar);

      filename = "/home/zouhair/left-or-right/a_" * string(α) * "_b_" * string(β) * ".jld"
      @printf("Saving %s: ", filename)
      @time save(filename, "Pistar", Pistar, "Ustar", Ustar, "α", α, "β", β);
    end
  end
end
@printf("Total run time for all cases: %.4f minutes", toc()/60)

end
