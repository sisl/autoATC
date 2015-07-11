module auxFuns

export isin


function isin{T}(x::T, L::Vector{T})
    for k in 1:length(L)
        if L[k] == x
            return true
        end
    end
    return false
end

end