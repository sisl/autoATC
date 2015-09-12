# README #

This repo contains Julia and Python code used as part of the autonomous air-traffic control for non-towered research.

Academic references available in the following two papers:

* Z. Mahboubi and M. J. Kochenderfer, “Autonomous Air Traffic Control for Non-Towered Airports,” in Air Traffic Management Research and Development Seminar, 2015. 
* Z. Mahboubi and M. J. Kochenderfer, “Continuous Time Autonomous Air Traffic Control for Non-Towered Airports,” in IEEE Conference on Decision and Control, Osaka, Japan, 2015. 

### How do I get set up? ###

* Install Julia (see http://julialang.org/downloads/) and download packages using Pkg.add(<PackageName>)
* Alternatively, run Julia remotely on JuliaBox (see https://juliabox.org/)

The project is organized as follows:

* model/ : contains files modelling the 3D dynamics of aircraft, and pilot behavior in airport pattern.
* model/solvers/   : different solvers used to find the optimal policies
* model/notebooks/ : iJulia notebooks used for development, visualization, post-processing, etc.
* MCTS_viz/ : used for visualizing trees resulting from Monte Carlo search
* modelLearning/ : different attempts at exploring empirical data for model learning (unpublished work)



### Who do I talk to? ###

* <zouhair@stanford.edu>