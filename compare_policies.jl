using CTMDP_kronsolver

loadPolicy(1.0f0, 0.0f0);

using CTMDP_mcts

function printStats(gammas, cnt, N)
    for γ in gammas
        println("γ=$γ -> state with different policies $(round(cnt[γ]/N*100,1))% [$(cnt[γ])]")
    end
    println("-------------------")
end



function sameAction(S::pattern.SType, p1::pattern.extActType, p2::pattern.extActType)
    if p1 == p2
        return true
    end
    
    ac1  = p1[1]; ac2  = p2[1];
    des1 = p1[2]; des2 = p2[2];
    
    if ac1 != 0 && ac2 != 0
        #Doesn't matter, we are just talking to different aircraft
        if S[ac1] == S[ac2] && des1 == des2
            return true
        end
    end
    
    return false
end



gammas = Float64[0.5, 0.7, 1.0]
cnt = Dict{Float64, Int64}()
for g in gammas cnt[g] = 0 end

dtprint = 60
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
        if !sameAction(S, p1, p2)
            #try again just in case since this is MCTS!
            p1 = mctsPolicy(S);
            if !sameAction(S, p1, p2)
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

