using CTMDP_kronsolver

loadPolicy(1.0f0, 0.0f0);

using CTMDP_mcts

function printStats(gammas, cnt, N)
    for γ in gammas
        println("γ=$γ -> state with different policies $(round(cnt[γ]/N*100,1))% [$(cnt[γ])]")
    end
    println("-------------------")
end

gammas = Float64[0.5, 0.7, 1.0]
cnt = Dict{Float64, Int64}()
for g in gammas cnt[g] = 0 end

dtprint = 6
start = time()
nextprint = start+dtprint
tic()

nVisited = 0
for cidx in 1:CTMDP_kronsolver.g_nScomp #rand(1:CTMDP_kronsolver.g_nScomp, N)
    S  = CTMDP_kronsolver.CIDX2S(cidx)
    p2 = ctmdpPolicy(S);
    nVisited += 1

    for γ in gammas
        mcts.pars.γ = γ
        p1 = mctsPolicy(S);
        if p1 != p2
            #try again just in case since this is MCTS!
            p1 = mctsPolicy(S);
                if p1 != p2 && p1[1] != 0 && !(string(S[p1[1]]) in ["LArr", "RArr"])
                cnt[γ] += 1
                #println(S, "-> M:", p1, " | C", p2)
                end
        end
    end
    now = time()
    if now > nextprint
        nextprint = now + dtprint
        printStats(gammas, cnt, nVisited)
    end
end
toc()
println("Final Stats:")
printStats(gammas, cnt, nVisited)

