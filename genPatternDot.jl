using airplaneType
using posType
using auxFuns
#############################################
##Printing to dot file for vizualization
#############################################

function drawNode(s, ned_start, ned_end)
    scale = 1000;
    w = (0.25*scale) / scale;
    off = (0.45*scale) / scale;
    
    x_s = ned_start.e / scale;
    y_s = ned_start.n / scale;
    x_e = ned_end.e / scale;
    y_e = ned_end.n / scale;
    
    
    str = """
    \\coordinate ($s S) at ($x_s,$y_s);
    \\coordinate ($s E) at ($x_e,$y_e);
    \\draw ($s S) -- ($s E);
    %\\path        ($s S) -- ($s E)
    %[
    %        sloped picture={
    %                \\draw [black] ($x_s,$y_s) rectangle ($x_e,$y_e) node[pos=0.5] {$s};
    %        }
    %];
    """
    
    
    angle = rad2deg(bearing(ned_start, ned_end))-90
    angle = rad2deg(atan2((y_e-y_s),(x_e-x_s)))
    dist = distance(ned_start, ned_end)/scale - off
    
    str = """
    \\draw[rounded corners,black,rotate around={$angle:($x_s,$y_s)}] ($(x_s),$(y_s-w)) rectangle ($(x_s+dist),$(y_s+w)) node[pos=0.5] {$s};
    
    """
    

   
   return str
        
end
        
function genPatternDot(posNE)
    fout = open("patternImg/patternImg.tex","w");
    
header = """
\\documentclass[tikz,border=2pt]{standalone}
\\makeatletter
\\usetikzlibrary{decorations.markings}
\\makeatother
\\begin{document}
\\tikzset{
        sloped picture/.style={
                decoration={
                        markings,
                        mark=at position 0.5 with {#1}
                },
                decorate
        }
}

\\begin{tikzpicture}
%\\draw [help lines] grid (3,2);
"""

    write(fout, header)
   
   
   
        
    for (s, ned_start) in posNE
        if s[2] == "E"
            continue
        end
        
        s_end = (s[1], "E")
        ned_end = posNE[s_end]
        str = drawNode(s[1], ned_start, ned_end)
        
        write(fout, str)  
    end
    
      
    
    
    
footer = """
\\end{tikzpicture}
\\end{document}
"""    
    write(fout, footer)
    close(fout)

    run(`pdflatex -halt-on-error -output-directory patternImg/ patternImg/patternImg.tex`)
end

#genPatternDot(airplaneType.posNE)



#############################################
##Printing to do file for vizualization
#############################################
function genPatternDotSimple(xy)
    fout = open("pattern.dot","w");
header="""
digraph G {
rankdir=LR

splines=true
node [shape = box]

"""
    write(fout, header)
    for f in keys(xy)
        sString = isin(f, pattern.sSafe) ? ", shape = doubleoctagon" : ""
        
        label = string(f)
        if f == :R
            label = "Runway (R)"
            sString *= ", width = 2"
        elseif f == :T
            label = "Taxi (T)"
            sString *= ", width = 1.5"
        end
        
    
        @printf(fout, "%s[label=\"%s\", pin=true, pos=\"%.2f, %.2f!\" %s]\n", f, label, xy[f][1], xy[f][2], sString)
    end
    
    for f in keys(xy)
        for t in pattern.NextStates[f]
            sString = ""
            ts = pattern.phaseFree(t)
            if f == :U2 && isin(ts, [:LDep, :RDep])
                sString = "[tailport=\"e\"]"
            end
            
            @printf(fout, "%s -> %s %s\n", f, ts, sString)
        end
    end
    
footer = "}"    
    write(fout, footer)
    close(fout)
end
genPatternDotSimple(pattern.xy)
run(`dot -Kfdp pattern.dot -Tpdf -o pattern.pdf`)