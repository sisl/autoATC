module sim3d_plotting

export plot3dsim, animate3dsim
export plotPerfResults


using airplaneType

using HDF5, JLD
using StatsBase

using PyPlot
using PyCall
@pyimport matplotlib.animation as anim
linewidth = 0.4
plt.rc("axes", linewidth=linewidth)
plt.rc("font", family="")
plt.rc("axes", titlesize="medium", labelsize="medium")
plt.rc("xtick", labelsize="small")
plt.rc("xtick.major", width=linewidth/2)
plt.rc("ytick", labelsize="small")
plt.rc("ytick.major", width=linewidth/2)
plt.rc("legend", fontsize="large")
plt.rc("axes" ,unicode_minus=false)


function plot3dsim(aircraftList, i, dL=20, dS=20; make3d = false)
    outOfData = false;
    b=0
    g=0
    if(make3d)
        ax = PyPlot.gca(projection="3d")
        #ax[:view_init](0, 180)
        ax[:pbaspect] = [1,1,0.6]   
    else
        #PyPlot.clf()
        ax = plt.gca()
            clf()
            hold("on")
    end
    
    D3 = Vector{Float64}[]
    XY = Vector{Float64}[]
        
    for ac in aircraftList
            range = (1:dL) + (i-1)*dS
            if(range[end] > length(ac.path))
                range = (length(ac.path)-dL):length(ac.path)
                outOfData = true 
            end
    
        Y = [p.n for p in ac.path[range]];
        X = [p.e for p in ac.path[range]];
        Z = -[p.d for p in ac.path[range]];
        N = length(X)
        if(make3d)
            di = 100
            for i in 1:di:(N-di)
                idx = i:(i+di)
                r = clip( mean(Z[idx]) / 300 , .05, 1.)    
                PyPlot.plot3D(Y[idx],X[idx],Z[idx],color=(r,g,0))
            end
            push!(D3, [Y[end], X[end], Z[end]])
            PyPlot.xlabel("North")
            PyPlot.ylabel("East")
            PyPlot.zlabel("Altitude")
        else
            PyPlot.plot(X,Y)
            #PyPlot.axis("equal")
            #PyPlot.ylabel("North")
            #PyPlot.xlabel("East")
        end
        g += 0.25
    end
    
        
    for (k, p) in airplaneType.posNE
        if k[2] != "S"
            continue
        end
        
        kk = (k[1], "E")
        
        pp = airplaneType.posNE[kk]
        
        x = [p.e, pp.e]
        y = [p.n, pp.n]
        
        plot(x,y,alpha=.5)
        text(mean(x), mean(y), string(k[1]),alpha=.3)
        hold("on")
    end
        
        
    for d in D3
        PyPlot.scatter3D(d[1],d[2],d[3],marker="o")         
    end
    #PyPlot.grid("on")
    if make3d
        ax[:set_ylim3d]([-5000, 5000])  
    else
        ylim([-1, 1]*6000)
        xlim([-1, 1]*12000)
        grid("on")
            hold("off")
    end
    
        
    return outOfData
    
end

function animate3dsim(aircraftList; nstart = 0, nend = Inf)

    # plot in an external window since it doesn't work yet in IJulia
    pygui(true)
    # tell PyPlot that the plot is interactive
    PyPlot.ion()
    # . . . and that previous plots are overwritten
    PyPlot.hold(false)
    # start time-stepping loop
    n = nstart;
    while(n < nend)
        outOfData = plot3dsim(aircraftList, n)
        # . . .
        # Then force the draw
        PyPlot.draw()
        n +=1 ;
        if(outOfData)
            break
        end
    end
    pygui(false)

end


function plotPerfResults(filelist)
    data = {JLD.load(f) for f in filelist}
    #TODO: find a better way to cycle through colors
    colors = ["b", "y", "c", "r", "g", "b", "k", "m","b"];
    for n in 1:length(data)
        atcIdx = 1
        
        nBetas = size(data[n]["flightTimes"],1)
        silentATCrun = size(data[n]["flightTimes"],2) == 2
        N = size(data[n]["flightTimes"],3)
    
        nmac_stats = zeros(nBetas, 2) #mean and std
        alert_stats = zeros(nmac_stats);
    
        ismcts = data[n]["policy"] == "MCTS"
        
        for b in 1:(nBetas)
            betaIdx = b
            ft = data[n]["flightTimes"][betaIdx,atcIdx,:][:]/3600
            nmacsPerHour  = ft ./ data[n]["nNMACcounts"][betaIdx,atcIdx,:][:]
    
            alertsPerHour = data[n]["alertCounts"][betaIdx,atcIdx,:][:] ./ ft
            #plt.plot(alertsPerHour/60, nmacsPerHour, linestyle="", marker=".")        
    
            nmacsPerHour = nmacsPerHour[!isinf(nmacsPerHour)] #Just in case we had things without incidents...
            nmac_stats[b,:]  = [mean_and_std(nmacsPerHour)...]
            alert_stats[b,:] = [mean_and_std(alertsPerHour/60)...]
            
        end
        
        xerrL = alert_stats[:,2]/sqrt(N)
        xerrR = alert_stats[:,2]/sqrt(N)
        yerrL = nmac_stats[:,2]/sqrt(N)
        yerrR = nmac_stats[:,2]/sqrt(N)
    
        
        prefix = ""
        pfix = "ctmdp"
        xpos = alert_stats[1, 1]+.05
        if ismcts
            prefix= ""
            pfix = "mcts"
            xpos -= .4
        end
        
        c = colors[n]
        
        plt.text(xpos , nmac_stats[1, 1]-.01,  prefix*string(data[n]["nPhases"])*"-phase-"*pfix, color=c)
    
        PyPlot.errorbar(alert_stats[:, 1], nmac_stats[:, 1], 
        xerr = [xerrL xerrR ]',  marker="o", linestyle = "--", yerr = [yerrL yerrR]', color=c)
        
    
    end
    PyPlot.xlim([-0.05, 1.8])
    PyPlot.grid("on")
    PyPlot.xlabel("Alert Rate (commands/minute)")
    PyPlot.ylabel("Mean Time to NMAC (hours)")
end

end