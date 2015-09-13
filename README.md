# README #

This repo contains Julia and Python code used as part of the autonomous air-traffic control for non-towered research.

Academic references available in the following two papers:

* Z. Mahboubi and M. J. Kochenderfer, “Autonomous Air Traffic Control for Non-Towered Airports,” in Air Traffic Management Research and Development Seminar, 2015. 
* Z. Mahboubi and M. J. Kochenderfer, “Continuous Time Autonomous Air Traffic Control for Non-Towered Airports,” in IEEE Conference on Decision and Control, 2015. 

### How do I get set up? ###

* Install Julia (see http://julialang.org/downloads/) and download packages using Pkg.add(<PackageName>)
    * Note that the code has only been tested with Julia Version 0.3.7
* Setup load paths by putting the following in ~/.juliarc.jl (if you don't want to dirty your LOAD_PATH, you can run those lines everytime julia is executed)
```
autoATCdir = "/path/to/autoATC/" 

push!(LOAD_PATH, autoATCdir * "model")
push!(LOAD_PATH, autoATCdir * "model/kron")
push!(LOAD_PATH, autoATCdir * "model/solvers")
```

### Project Structure ###

The project is organized as follows:

* model/ : contains files modelling the 3D dynamics of aircraft, and pilot behavior in airport pattern.
    *  solvers/   : different solvers used to find the optimal policies
    *  notebooks/ : iJulia notebooks used for development, visualization, post-processing, etc.
    *  patternImg/ : used to generate the locations of the pattern in the model
* modelLearning/ : different attempts at exploring empirical data for model learning (unpublished work)
* MCTS_viz/ : used for visualizing trees resulting from Monte Carlo search


Additional notes:

The number of phases (nPhases) used in the model for the pattern is defined in `model/pattern.jl`. Unfortunately, 
this needs to be modified on a case-by-case as Julia does not make it easy to define parameters in modules
at runtime. As a workaround, the file patternSweep.sh can be used to run the same script using multiple
values for nPhases.

The file `model/solvers/solveCTMDPs.jl` solves the continous-time MDP formulation. It relies on the CTMDP_kronsolver module. 

The file `parallelRun.jl` (wrapper for `runSims_paralle.jl`) is used to simulate different policies using processes in parallel.
Note that depending on policy type, policies might need to be solved ahead of time.


### Who do I talk to? ###
The code does not come with any support or guarantees. However, feel free to contact the author if you have questions or suggestions.