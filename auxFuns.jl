module auxFuns

export isin, swap, swap!


function isin{T}(x::T, L::Vector{T})
    for k in 1:length(L)
        if L[k] == x
            return true
        end
    end
    return false
end





function swap{T}(s::Vector{T}, i, j)
  sc = copy(s);
  if(i != 0 && j != 0)
    sc[i] = s[j];
    sc[j] = s[i];
  end

  return sc
end

function swap!{T}(s::Vector{T}, i, j)
  if(i != 0 && j != 0)
    tmp = s[i]
    s[i] = s[j];
    s[j] = tmp;
  end
end




end


