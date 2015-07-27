module CTMDP_kronindexing

using pattern
using auxFuns

export CIDX2S, CIDX2X, X2CIDX, S2CIDX
export LIDX2CIDX, CIDX2LIDX

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
# No tail recursion optimization in julia :(
# function largest!(i::Int64, n::Int64, k::Int64, val_off::Vector{Int64})
#    x = binomial(n+k-1, k)   #g(n,k) 
#    if x <= i
#         @inbounds val_off[1] = n
#         @inbounds val_off[2] = x
#         return nothing
#    end
#    #tail recursion!
#    return largest!(i, n-1, k, val_off)
# end

function largest!(i::Int64, n::Int64, k::Int64, val_off::Vector{Int64})
   x = binomial(n+k-1, k)
   
   while(x > i)
       n -= 1
       x = div(x*n, n+k) #avoids recalling binomial(n+k-1,k) !!
   end
   
   @inbounds val_off[1] = n
   @inbounds val_off[2] = x
   
   return nothing
end

#Consider reusing this and rely on the fact that nothing is multithreaded!???
const val_off = Int64[0,0]

function id2combo!(combo::Vector{Int64}, idx::Int64, n::Int64, k::Int64)
    #This is based on http://stackoverflow.com/questions/12146910/finding-the-lexicographic-index-of-a-permutation-of-a-given-array
    
    #We assume combo that idx is 1 based indexing
    #transform it to 0 based indexing
    idx -= 1 #Ooo the joys of 1 based indexing ...

    #val_off = Int64[0,0] 
    while(k >= 1)
        largest!(idx,n,k,val_off)
        @inbounds combo[k] = val_off[1] + 1 #Ooo the joys of 1 based indexing!
        @inbounds offset = val_off[2]
        k -= 1
        idx -= offset
    end
    return nothing
end

function combo2id(combo)
    #This assumes combo is 1 based and sorted in increasing order
    
    #see http://stackoverflow.com/questions/18613690/calculate-nth-multiset-combination-with-repetition-based-only-on-index
    #the ϕ bijection goes from combo[i] -> combo[i]+i
    #the f bijection is given on wiki http://en.wikipedia.org/wiki/Combinatorial_number_system
    #All of this assumed 0 based indexing. So we need to be careful with Julia's 1 based indexing ...
    
    id = 0 #0 based index 
    k = length(combo)
    for i1 in 1:k #1 based indexing
        i0 = i1 - 1 #0 based indexing
        id += binomial((combo[i1]-1)+i0, i1)
    end
    return id+1 #convert back to 1 indexing
end
###########################################
###########################################

#Going from compact indices to states
function CIDX2S(cidx::Int64)
  return X2S(CIDX2X(cidx))
end
function CIDX2X(cidx::Int64)
  combo = zeros(xType, g_nVehicles)
  return CIDX2X!(combo, cidx)
end
function CIDX2X!(combo::XType, cidx::Int64)
  id2combo!(combo, cidx, g_nNodes , g_nVehicles)
  return combo
end

#Going from states to compact indices
function X2CIDX(X::XType)
  return combo2id(sort(X))
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
  X = CIDX2X(cidx)
  return X2LIDX(X) 
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

end