if isdefined(current_module(),:LastMain) && isdefined(LastMain,:Iterators)
  using LastMain.Iterators
else
  using Iterators
end


#Number of instances
g_nI = 2
#Number of nodes per instaces
g_N = 5


function combos_with_replacement(list, k)
    n = length(list)
    [[list[c[i]-i+1] for i=1:length(c)] for c in combinations([1:(n+k-1)],k)]
end

#Get all possible states, allowing replacement, but order does not matter:
const g_S = combos_with_replacement(1:g_N, g_nI)
const g_NS = length(g_S)
#To avoid having to hash this all the time, we'll create an S2i function
const g_S2i_dict = Dict(g_S,[1:g_NS])
function S2i(s::typeof(g_S[1]))
  return g_S2i_dict[sort(s)]
end


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


function swap(s, i, j)
  sc = copy(s);
  if(i != 0 && j != 0)
    sc[i] = s[j];
    sc[j] = s[i];
  end

  return sc
end

function s2Qidx(s::Vector{Int64},N::Integer)
  Ns = length(s)
  return sub2ind(N*ones(Int64, Ns), reverse(s)...)
end
function QVeval(s, action, Qdict,V)
  #We'll assume that the actions
  #are (idx, act), Qdict is
  (idx, act) = action;
  #First, get the Q that is relevant
  actQ = [:∅ for i in 1:g_nI]
  actQ[1] = act;
  Q = Qdict[actQ]

  #Next, re-order the
  #states so that we are addressing
  #the right instance

  so  = swap(s, 1, idx)
  sop = swap(sp, 1, idx)


  #The last thing left to do is figure out
  #how to index into Q based on the ordering!
  si  = s2Qidx(so, g_N)
  spi = s2Qidx(sop, g_N)

  #Finally, we are ready to evaluate Q and return
  #the value
  return Q[si, spi]
end
