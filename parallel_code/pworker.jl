# parallel test code

global __PARALLEL__

if !isdefined(:__PARALLEL__)
    __PARALLEL__ = false
end

if !__PARALLEL__
    # do not use PyPlot in parallel
    using PyPlot
end

using Distributions


function runTest(n::Int64; bParallel::Bool = false)

    ID = zeros(Uint64, n)
    N = zeros(Uint64, n)

    u = 200000
    lsts = rand((u-int64(u*0.1)):(u+int64(u*0.1)), n)

    if !bParallel
        k = 1

        for x = lsts
            result = countPrimes(x)

            ID[k] = result[1]
            N[k] = result[2]

            k += 1
        end

    else
        results = pmap(x -> countPrimes(x), lsts)

        k = 1
        for result in results
            ID[k] = result[1]
            N[k] = result[2]
            k += 1
        end
    end

    for i = 1:n
        println("id = ", ID[i], ", x = ", lsts[i], ", n_primes = ", N[i])
    end
end


function countPrimes(x::Int64)

    nPrimes = 0
    n = 0

    for i = 2:x
        bPrime = true

        if i > 2
            for j = 2:i-1
                n += 1

                if i % j == 0
                    bPrime = false
                    break
                end
            end

        end

        if bPrime
            #println(i)
            nPrimes += 1
        end
    end

    return myid(), nPrimes
end


if !__PARALLEL__
    runTest(4)
end


