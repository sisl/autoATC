module CTMDP_kron

using pattern
using auxFuns
using CTMDP_kronindexing

using HDF5, JLD


export loadCTMDPpolicy, saveCTMDPpolicy
export ctmdpPolicy


function policy_X2a_compact(X::XType, Aopt::Vector{compActType})
    Xperm = sortperm(X)
    X_cidx = X2CIDX(X)
    
    compactAct = Aopt[X_cidx]
    
    act = pattern.g_nullAct
    if compactAct != pattern.g_nullAct
      pidx = int8(Xperm[compactAct[1]])
      act = [pidx, compactAct[2]]
    end
    
    #This is still in compact form
    return act
end
##################################
const Aopt = Array(compActType, g_nScomp);
for i in 1:g_nScomp
    Aopt[i] = copy(g_nullAct)
end

function ctmdpPolicy(S::SType)
    X = S2X(S)
    #Get the compact form representation, accounting for permutation
    act = policy_X2a_compact(X, Aopt::Vector{compActType})
    #Trasnform it to the extended form
    return pattern.compAct2extAct(act,S)
end

##################################
function saveCTMDPpolicy(Aopt, α, β_cost; prefix="")
    filename = "policies/" * prefix * "CTMDPpolicy_n_" * string(pattern.nPhases) * "_a_" * string(α) * "_b_" * string(β_cost) * ".jld"
    Aopt_idx = Array(Int8, g_nScomp)
    Aopt_act = Array(Int8, g_nScomp)
    for i in 1:g_nScomp
        Aopt_idx[i] = Aopt[i][1]
        Aopt_act[i] = Aopt[i][2]
    end 
    JLD.save(filename, "Aopt_idx", Aopt_idx, "Aopt_act", Aopt_act); #"Vstar", Vshort,
    println(filename)
end

function loadCTMDPpolicy(α, β_cost; prefix="")
    filename = "policies/" * prefix * "CTMDPpolicy_n_" * string(pattern.nPhases) * "_a_" * string(α) * "_b_" * string(β_cost) * ".jld"
    data = JLD.load(filename)
    global Aopt

    assert(length(data["Aopt_idx"]) == length(Aopt))
    #These are just references
    Aopt_idx = data["Aopt_idx"]
    Aopt_act = data["Aopt_act"]
    for i in 1:g_nScomp
        Aopt[i][1] = Aopt_idx[i]
        Aopt[i][2] = Aopt_act[i]
    end 
    return ctmdpPolicy
end

end