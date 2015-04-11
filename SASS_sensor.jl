module SASS_sensor

using posType
using pattern
using airplaneType

export ne_psi, SASS_sense!

#Convert sigma=10ft and .4 degs to SI units
const σ_ρ = 10 / 3.33 #10 feet according to SASS report from MIT
const σ_θ = .4 * pi/180. # 0.4 degs according to SASS report from MIT

#Noise in heading measurement
const σ_ψ = 10*σ_θ


const baseLine = 2000/3
const sensors_pos = [pos(0, baseLine/10, 0.) , pos(baseLine, -baseLine/10, 0.)]
const nSensors = length(sensors_pos)
const frac = 1./nSensors;
#Measurement model
#This function returns x,y position of the aircraft
#by determining the 


function theta(p0::pos, p1::pos)
    dN = p1.n - p0.n;
    dE = p1.e - p0.e;
    return atan2(dN, dE)
end


type ne_psi
    n::Float64
    e::Float64
    psi::Float64
    function ne_psi(n,e,psi)
        new(n,e,psi)
    end
    ne_psi() = ne_psi(0.,0.,0.)
end

function SASS_sense!(p_hat::ne_psi, ac::airplane)
    p_hat.e = 0.
    p_hat.n = 0.
    for i in 1:nSensors
        #compute bearing and range to each sensor and add noise
        θ_hat = theta(sensors_pos[i], ac.posNED)     + σ_θ*randn(pattern.rng)
        ρ_hat = distance(sensors_pos[i], ac.posNED)  + σ_ρ*randn(pattern.rng)
        #Project back to get position. Average all sensors
        p_hat.e += (sensors_pos[i].e + ρ_hat * cos(θ_hat))*frac
        p_hat.n += (sensors_pos[i].n + ρ_hat * sin(θ_hat))*frac
    end
    
    #return the heading of the aircraft with appropriate noise
    p_hat.psi = unwrap(ac.psi + σ_ψ*randn(pattern.rng))
    return nothing
end

end