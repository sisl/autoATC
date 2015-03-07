include("pattern.jl")
include("CTMDP.jl")
β_cost = 0.0;
legalActions = validCompactActions
amax = maximum([length( NextStates[s]) for s in keys(NextStates)]);
A = (Int64)[0:amax];
P0 = spzeros(g_nNodes, g_nNodes);
P = (typeof(P0))[copy(P0) for a in A]

for a in A
    if a == 0
        act = g_nullAct
    else
        act = Int64[1, a]
    end
    for x in 1:g_nNodes
        s = x2s(x)
        for sp in NextStates[s]
            xp = s2x(sp)
            P[a+1][x,xp] = Transition([s], act, [sp])
        end
    end
end
M0 = speye(g_nNodes)
for x in 1:g_nNodes
    s = x2s(x)
    M0[x,x] = 1./(teaTime[s]/60)
end

Isp = spdiagm(ones(g_nNodes));
Qt_list = (typeof(M0))[(M0*(P[a+1] - Isp))' for a in A];
Vshort = zeros(g_nScomp);

#Aopt = gaussSeidel!(Qt_list, Vshort, 0.1);# at /home/zouhair/Dropbox/Research/autoATC/CTMDP.jl:268
