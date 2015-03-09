###############################
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

###############################
#Kronecker vectors and matrices
function ev(i, m)
    e = spzeros(m,1)
    e[i] = 1.;
    return e;
end

function E(i,j,m,n)
    #Same as Em = ev(i,m) * ev(j,n)'
    Emn = spzeros(m,n)
    Emn[i,j] = 1.
    return Emn
end

###############################
function Ppq_v!(p,q, v)
    assert(issparse(v) && size(v,2)==1)
    #equivalent to v = Pcrazy(p,q) * v
    for i in v.colptr[1]:(v.colptr[2]-1)
        l = v.rowval[i]
        #k = sub2ind((p,q), reverse(ind2sub((q,p),l))...)
        #Above allocates memory and needs to reverse. This one liner is faster
        #v.rowval[i] = rem(l-1,q) * p + div(l-1,q) + 1
        #The one below is even faster (4x compared to original!) 
        (d,r) = divrem(l-1,q)
        v.rowval[i] =  r * p + d + 1
    end
    return v
end

function Pmn(m,n)
    Pmn_ret = spzeros(m*n,m*n)
    for i in 1:m
        for j in 1:n
            k = sub2ind((n,m),j,i)
            l = sub2ind((m,n),i,j)
            Pmn_ret[k,l] = 1.
        end
    end
    return Pmn_ret
    #Slow version
    #Pmn = spzeros(m*n,m*n)
    #for i in 1:m
    #    for j in 1:n
    #        Eij = E(i,j,m,n)
    #        Pmn = Pmn + kron(Eij, Eij')
    #    end
    #end
    #return Pmn
end

###############################
#kron(ea,v)
###############################
# function eaKronv(a,n,v)
#     m = length(v)
#     res = spzeros(n*m,1)
#     idx = m*(a-1) + (1:m)
#     res[idx] = v
#     return res
# end
#In place version.
function eaKronv!(a,n,v)
    offset = v.m*(a-1)
    v.m = v.m * n;
    for i in v.colptr[1]:(v.colptr[2]-1)
        v.rowval[i] =  offset + v.rowval[i]
    end
end

###############################
#kron(v,ea)
###############################
# function vKronea(v,a,n)
#     m = length(v)
#     res = spzeros(n*m,1)
#     idx =  a:n:(n*m)
#     res[idx] = v
#     return res
# end
#In place version.
function vKronea!(v,a,n)
    v.m = v.m * n;
    for i in v.colptr[1]:(v.colptr[2]-1)
        v.rowval[i] =  a + (v.rowval[i]-1)*n
    end
end


function sparseVectorPlusEq!(u, v, n)
    #Add the B values in to A
    for i in v.colptr[1]:(v.colptr[2]-1) 
        u[v.rowval[i]] += v.nzval[i]
    end
end

###############################
#This is the heart of most of it
function Cbt(Bt,K,b)
    n = size(Bt,1);
    n_K = n^K;
    n_Km1 = n^(K-1);

    #Compute Cb
    Cb_res = spzeros(n_K,1);

    #Note that we will abuse this vector!
    res_u = spzeros(n, 1)
    #n should be a manageable size...
    #might consider moving this out of this function though?
    res_u_rowval = Array(Int64, n)
    res_u_nzval = Array(Float64, n)
    for u in 0:(K-1)
        p = n^u;
        q = n^(K-u)
        bp = sub2ind((p,q), reverse(ind2sub((q,p),b))...)
        (d,c) = ind2sub((n,n_Km1), bp)
        
        #reisze res_u and populate it
        res_u.m = n
        #This for loop grabs the d column of Bt
        idx = 1
        for i in Bt.colptr[d]:(Bt.colptr[d+1]-1) 
          res_u_rowval[idx] = Bt.rowval[i]
          res_u_nzval[idx] = Bt.nzval[i]
          idx += 1
        end
        res_u.colptr[2] = idx
        res_u.rowval = res_u_rowval
        res_u.nzval = res_u_nzval
        eaKronv!(c, n_Km1, res_u)
        #Last step is to permute
        Ppq_v!(q,p,res_u)
        #Cumulative sum into accumulator
        #Cb_res = Cb_res + res_u
        sparseVectorPlusEq!(Cb_res, res_u, n_K)
    end

    return Cb_res
end
###############################
#This is the function that puts it all together!

function Qti(A,B,K,i)
    return Qi_ABt(A',B',K,i)
end
function Qti_ABt(At,Bt,K,i)
    n = size(At,1)
    assert(n == size(At,2)) #enforce squareness
    assert(size(At) == size(Bt)) #only working with same size matrices

    n_K = n^K;
    (b, a) = ind2sub((n_K , n), i);

    resA = At[:,a];
    vKronea!(resA,b,n_K)

    res = Cbt(Bt,K,b)
    eaKronv!(a,n,res)
    
    res = res + resA
    #sparseVectorPlusEq!(res, resA, n_K)


    return res
end

#This is the lazy version
#pretty slow, and will run out
#out of memory for large n's
function Qti_lazy(A,B,K,i)
    Q = A;
    for k in 1:K
       Q = kronSum(Q,B)
    end

    return Q[i,:]'
end
