# TRR-ACoP2022-Public
Public release version of ACoP2022

This is the code for our ACoP-2022 poster:
An Integrated QSP Model of Lipoproteins with Liver Metabolism for NAFLD
Authors: Theodore Rieger, Anna Sher, Tyler Cassidy, and C.J. Musante

The code was written and tested in Julia v1.7.3 and v1.8.5

Execution: `include("start.jl")`

Note: start.jl will add and pre-compile all the dependencies, thus execution will be very slow the first time.

The final poster and target figures are included in the repository.

Folder structure:
* fig/ -- Figure files, pre-populated with the expected output from the scripts.
* lit/ -- Useful literature
* log/ -- Initial parameter estimate workbooks and PPop jld2 files, cleanup on this directory is a TODO
* public/ -- Public data from NHANES
* src/ -- Main source code directory, described below

src/ contents:
* drawp.jl -- Draws a new set of parameters from a passed distribution.
* errorellipse.jl -- Returns a 2D error ellipse at the 95% prediction interval, used for plotting.
* figure1to3.jl -- Long script for actually plotting the final figures. Note: despite the name, this plots 4 figures. :)
* mh.jl -- Home-brewed Metropolis-Hastings algorithm that is friendly to QSP models. Based on Rieger et al. 2018.
* mh_sim.jl -- Helper script for MH algorithm, all problem specific information is here, including the objective function.
* model.jl -- Main script, called by start.jl. Currently, the script uses the pre-built PP files in the log/ directory for simulation and analysis.
* readnhanes.jl -- Helper script for ingesting and fitting the NHANES dataset.
* rxns.jl -- Reaction definitions of the model.
* select_vps.jl -- Acceptance-rejection sampling of Allen et al. 2016. Translated from Matlab.
* therapy_ebs.jl -- Calculation of the error-bounds for the literature data (percent changes).


