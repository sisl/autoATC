
importall posType

#Convert sigma=10ft and .4 degs to SI units
const σ_ρ = 10 / 3.33 #10 feet according to SASS report from MIT
const σ_θ = .4 * pi/180. # 0.4 degs according to SASS report from MIT


const baseLine = 2000/3
const sensor1_pos = pos(0, baseLine/10, 0.)
const sensor2_pos = pos(baseLine, -baseLine/10, 0.)

#Measurement model
#This function returns x,y position of the aircraft
#by determining the 

