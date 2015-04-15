module posType

export pos, bearing, distance2, distance, project, unwrap, clip

#################################################
type pos
#################################################
  n::Float64
  e::Float64
  d::Float64
  
  function pos(n,e,d)
    new(n,e,d)
  end
  pos() = pos(0., 0., 0.)
  pos(ned::Vector{Float64}) = pos(ned[1], ned[2], ned[3])
end
#################################################

function bearing(p0::pos, p1::pos)
    dN = p1.n - p0.n;
    dE = p1.e - p0.e;
    return unwrap(pi/2-atan2(dN, dE))
end


function distance2(p0::pos, p1::pos)
    dN = p0.n - p1.n;
    dE = p0.e - p1.e;
    return (dN*dN + dE*dE)
end

function distance(p0::pos, p1::pos)
    return sqrt(distance2(p0,p1))
end


function project(p0::pos, distance, bearing)
  dN =  distance * cos(bearing)
  dE =  distance * sin(bearing)
  return pos(p0.n + dN, p0.e + dE, p0.d);
end


function +(p0::pos, p1::pos)
    return pos(p0.n + p1.n , p0.e + p1.e, p0.d + p1.d)
end


import Base.copy
function copy(p0::pos)
    return pos(p0.n, p0.e, p0.d)
end



function unwrap(a::Float64)
    if (a > pi)
        return mod(a + pi, 2*pi) - pi
    elseif (a < -pi)
        return -(mod(-a + pi, 2*pi) - pi)  #aka -unwrap(-a)
    else
        return a;
    end
end

function clip(r::Float64, min::Float64, max::Float64)
  if r < min
    return min
  elseif r > max
    return max
  else
    return r
  end
end




end