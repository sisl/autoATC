if isdefined(current_module(),:LastMain) && isdefined(LastMain,:Iterators)
  using LastMain.Iterators
else
  using Iterators
end


#Number of instances
g_nI = 2
#Number of nodes per instaces
g_N = 5

g_noaction = (0, :∅)


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


function findn_rows{Tv,Ti}(S::SparseMatrixCSC{Tv,Ti}, colIdx::Integer)
    idx = S.colptr[colIdx] : (S.colptr[colIdx+1]-1)
    return (S.rowval[idx] , S.nzval[idx])
end


function QVeval(s, action, Qtdict, V::Vector{Float64}, β::Float64)
  #We'll assume that the actions
  #are (idx, act), Qtdict is
  (idx, act) = action;
  #First, get the Q that is relevant
  actQ = (Symbol)[:∅ for i in 1:g_nI]
  actQ[1] = act;
  #TODO: Make sure this is a reference
  # and NOT a copy!
  Qt = Qtdict[actQ]

  #Next, re-order the
  #states so that we are addressing
  #the right instance
  so  = swap(s, 1, idx)


  #The last thing left to do is figure out
  #how to index into Q based on the ordering!
  si  = s2Qidx(so, g_N)

  #These are the rows of Qt for this column,
  #which means these are the columns of Q for this
  #state (s -> si). These indices can be used to
  #index into the value function V!
  indices = Qt.colptr[si] : (Qt.colptr[si+1]-1)
  vIndices = Qt.rowval[indices]

  #Diagonal element is negative by construction
  #But q(x) is defined as the positive value
  qx = -Qt[si,si]

  #q(x) is not supposed to be part of the summation
  #we will later add Qt[si,si]*V[si] so starting with -Qt[si,si] will cancel it out
  qVsum =  qx * V[si]

  for idx in indices
    Qval = Qt.nzval[idx]

    vIdx = Qt.rowval[idx]
    #TODO: take advantage of the state ordering for V!
    qVsum += Qval * V[vIdx]
  end

  return (qVsum + r(s, action)) / (β + qx)
end

#TODO: Implement this properly for
#autoATC
#Consider precomputing r ?
function r(s::Vector{Int64}, a)
  #rate in state 1 is the best
  #unless two things are in the same state!

  Ns = length(s)
  Ns_u = length(unique(s))

  R = 0.

  if(Ns != Ns_u)
    R = -10000.
  elseif 1 in s
    R = 100.
  end

  if(a != g_noaction)
    R -= 100.
  end

  return R;
end
