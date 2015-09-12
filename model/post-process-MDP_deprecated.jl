Pistar = Array(Int32, NS)


function randomStart()
  rstart = ()-> allstates[rand(3:length(allstates), length(S[1]))];
  newS0 = rstart()
  while(NcolNtaxi(newS0)[1] != 0)
    newS0 = rstart()
  end
  return newS0
end



α_range = [0., 0.05, 0.25, 0.5, 0.75, 0.85, 0.95, 1.]
β_range = [-1., α_range] #-1 is used for the alpha analysis to use the true alpha
#β_range = [0., 0.01, 0.1, 0.25, 0.5, 1.]

meanCollisions = zeros(Float32, length(α_range), length(β_range))
stdCollisions = zeros(Float32, length(α_range), length(β_range))
firstCollisions = zeros(Float32, length(α_range), length(β_range))
verb = zeros(Float32, length(α_range), length(β_range))


runAlphaAnalysis = true
srand(rng, uint32(3))
for i in 1:length(α_range)
  α = α_range[i];
  for j in 1:length(β_range)
    if runAlphaAnalysis
      β = 0.5
      fakeAlpha = β_range[j]
    else
      β = β_range[j];
    end

    #@printf("α = %.2f, β = %.2f\n", α, β)

    #fileName = "/home/zouhair/left-or-right/a_" * string(α) * "_b_" * string(β) * ".jld"
    if !(runAlphaAnalysis) || fakeAlpha == -1.
      fileName = "/home/zouhair/left-or-right/a_" * string(α) * "_b_" * string(β) * ".jld"
    else
      fileName = "/home/zouhair/left-or-right/a_" * string(fakeAlpha) * "_b_" * string(β) * ".jld"
    end
    d = load(fileName)

    Pistar = d["Pistar"]
    verb[i,j] = countnz(Pistar-1)/(1. * NS)

    #simulate 20x 1000's starting from random locations
    ncol = Array(Float32, 20);
    tcol = Array(Float32, 20);
    for n in 1:length(ncol)
      S0 = randomStart()
      (Satc, Aatc, natc, ratc) = simulate(S0, PiStarFun, N=1000);
      ncol[n] = cumsum(natc)[end] / 1000.
      if(ncol[n] > 0)
        tcol[n] = find(natc)[1] #time of first collision
      else
        tcol[n] = 1000.;
      end
    end
    meanCollisions[i,j] = mean(ncol);
    stdCollisions[i,j] = std(ncol)

    firstCollisions[i,j] = mean(tcol)

  end
end
