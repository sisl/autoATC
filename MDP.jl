
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
#Policy function
#############################################
function PiStarFun(s0::Array{Symbol, 1})
    s_i = s2i[s0];
    return validActs[s_i][Pistar[s_i]]
end

function PiStarFunSlow(s0::Array{Symbol, 1})
    s_i = s2i[s0];
    return validActions(s0)[Pistar[s_i]]
end



if(false)


@printf("Allocating memory for transitions/Rewards \n"); tic();
tsaEntry = typeof(Float32[])[]
rsaEntry = Float32[];
TSA = zeros(typeof(tsaEntry), NS)
RSA = zeros(typeof(rsaEntry), NS)
Jlist = zeros(typeof(Int64[]), NS)
@time for from_i in 1:NS
  if(from_i % int(NS/4) == 0)
    toc(); @printf("%i/%i ... \n",from_i, NS); tic();
  end

  from = S[from_i]
  nn =  length(validActs[from_i])

  RSA[from_i] = zeros(Float32, nn)
  TSA[from_i] = zeros(typeof(Float32[]), nn)

  from_i_nextPerms = NextPerms[from_i]
  Jlist[from_i] = Int64[s2i[to] for to in from_i_nextPerms]

  for (a_i, a) in enumerate(validActs[from_i])
      TSA[from_i][a_i] = zeros(Float32, length(from_i_nextPerms))
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
    Piarray = zeros(Int64, NS)

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
#Solving
#############################################
if !isdefined(current_module(),:solveAll)
  solveAll = false
end


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
