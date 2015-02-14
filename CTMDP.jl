if isdefined(current_module(),:LastMain) && isdefined(LastMain,:Iterators)
  using LastMain.Iterators
else
  using Iterators
end

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


function swap(s, i, j)
  sc = copy(s);
  if(i != 0 && j != 0)
    sc[i] = s[j];
    sc[j] = s[i];
  end

  return sc
end



function combos_with_replacement(list, k)
    n = length(list)
    [[list[c[i]-i+1] for i=1:length(c)] for c in combinations([1:(n+k-1)],k)]
end

###########################################
#problem definitions
###########################################


#Number of instances
g_nVehicles = 2
#Number of nodes per instaces
g_nNodes = 5

g_noaction = (0, :∅)



#Get all possible states, allowing replacement, but order does not matter:
const g_Sshort = combos_with_replacement(1:g_nNodes, g_nVehicles)
const g_nSshort = length(g_Sshort)
#To avoid having to hash this all the time, we'll create an s2shortidx function
const g_s2shortidx_dict = (typeof(g_Sshort[1]) => Int64)[g_Sshort[sidx] => sidx for sidx in 1:g_nSshort]


###########################################
#Awesome functions for indexing magic :)
###########################################

function s2shortidx(s::typeof(g_Sshort[1]))
  return g_s2shortidx_dict[sort(s)]
end
function shortidx2s(sidx::Int64)
  return g_Sshort[sidx]
end

function s2longidx(s::Vector{Int64})
  #Note the reverse due to the nature of s2longidx!
  return sub2ind(g_nNodes*ones(s), reverse(s)...)
end

function longidx2s(lidx::Int64)
  dims = g_nNodes*ones(Int64, g_nVehicles)
  s_reversed = [ind2sub(tuple(dims...), lidx) ...]
  return reverse(s_reversed)
end

function lidx2sidx(lidx::Int64)
  dims = g_nNodes*ones(Int64, g_nVehicles)
  s_reversed = [ind2sub(tuple(dims...), lidx) ...]
  #Note how we don't bother about reversing
  #things since s2shortidx sorts things anyway!
  return s2shortidx(s_reversed)
end

#Note that there is not a unique lidx
#that corresponds to a given sidx!
function sidx2lidx(sidx::Int64)
  return s2longidx(g_Sshort[sidx])
end
###########################################


function findn_rows{Tv,Ti}(S::SparseMatrixCSC{Tv,Ti}, colIdx::Integer)
    idx = S.colptr[colIdx] : (S.colptr[colIdx+1]-1)
    return (S.rowval[idx] , S.nzval[idx])
end


function QVeval(s, action, Qtdict, Vshort::Vector{Float64}, β::Float64)
  #We'll assume that the actions
  #are (idx, act), Qtdict is
  (idx, act) = action;
  #First, get the Q that is relevant
  actQ = (Symbol)[:∅ for i in 1:g_nVehicles]
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
  si  = s2longidx(so)

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
  qVsum =  qx * Vshort[s2shortidx(so)]

  for idx in indices
    Qval = Qt.nzval[idx]

    vIdx = Qt.rowval[idx]
    qVsum += Qval * Vshort[lidx2sidx(vIdx)]
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


function gaussSeidel!(Vshort)
  Aopt = [g_noaction for i in 1:g_nSshort];

@time for i in 1:100
    maxVchange = 0.
    for (si, s) in enumerate(g_Sshort)
        β = 1./0.95
        aopt = g_noaction
        Qmax = QVeval(s, g_noaction, Qt_joint, Vshort, β)
        for a in product([1,2], [:L, :R, :S]) #should be validActions(s)
            Qa = QVeval(s, a, Qt_joint, Vshort, β)
            if Qa > Qmax
                Qmax = Qa
                aopt = a
            end
        end

        maxVchange = max(maxVchange, abs(Vshort[si] - Qmax))
        Vshort[si] = Qmax
        Aopt[si] = aopt
    end

    if(maxVchange < 1)
        println("stopping after ", i)
        break
    end
end

return Aopt

end
