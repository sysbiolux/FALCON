FALCON
Fast algorithm for the contextualisation of logical networks.

Sébastien De Landtsheer, sebastien.delandtsheer@uni.lu

::: Contents :::

ExampleDatasets	: collection of examples from published papers (De Landtsheer et al., 2016, Bioinformatics, De Landtsheer et al., 2017) to reproduce and benchmark the toolbox

ThirdParty	: collection of scripts from other sources, used by some functions of the toolbox

Source	: the actual toolbox. Contains the following scripts:

FalconInstall.m     : add Falcon scripts and example to the path (type FalconInstall on the Command Window)

FalconMakeModel.m   : creates the optimization problem from a list of interactions (InputFile) and a list of measured nodes in different conditions (MeasFile). Network files can be txt (tab-delimited, old deprecated format) or .sif or .xls or .xlsx, while data files can be .txt (tab-delimited, old deprecated format) or .xls, or .xlsx, or .csv. See examples in 'PDGF' forlder for format examples. Automatically assigns the weights of interactions for single-input nodes and creates the constrains for the optimization problem.

FalconMakeGlobalModel : creates an optimization problem as above, for a number of different contexts. A meta-model is created for which each context will be modeled independently. Parameters can be regularized accross contexts, cfr infra.

FalconExpand        : automatic expansion of logical networks in order to facilitate rapid computation

FalconInt2File.m    : write the list of new interactions into a temporary file for re-importing (essential to generate perturbed model variants)

FalconIC.m          : automatically assigns initial parameter values which are drawn from a uniform or normal distribution.

FalconObjFun.m      : optimization with built-in objective function that simulates networks and returns the sum of squared errors between the simulated and expected values. Uses fmincon with the default algorithm (interior-point). Can be used inside a 'parfor' loop to perform multiple optimizations in parallel.

FalconObjFunPlus.m  : same as FalconObjFun.m but returns a global variable containing the costs over the course of the optimization. Useful for convergence diagnostics. Cannot be used inside a parfor loop.

FalconObjFun_mute 	: same as FalconObjFun.m but only returns a low amount of info to console during execution. 

FalconResults.m     : summarises the optimization output and display the best and average values for the parameters.

FalconFitEvol.m     : re-optimise the model and plot the evolution of fitting costs (3 rounds)

FalconSimul.m       : resimulates the optimized network,displays the results as plots, and saves a bunch of figures and tables to the result folder.

FalconShowNetwork.m : uses the 'biograph' environment to display the interaction graphs of the networks.
***known issue with this function on some platforms***

FalconLPSA.m        : performs a local parameters sensitivity analysis on an optimized model.

FalconKO.m          : creates new models for each knock-outed parameter by setting parameter value to 0 and re-optimise

FalconKONodes.m     : FalconKO creates new models for each knock-outed node by setting outcoming parameter value(s) to 0 and re-optimise

FalconResample.m    : creates a new dataset from the original one, given the mean and standard deviation of the original one. Useful to assess the variability of the parameter estimates, in conjunction with FalconLPSA.m. ***note: more options are available with the new function FalconPreProcess.m (cfr infra)

FalconPreProcess	: see in-function help for 

FalconValidation.m  : resimulates the optimized network and compare to the validation dataset

FalconPlots.m       : re-plot the figures from the chosen analysis (in case the plot(s) crashed during the analytical process)

FalconGUI.m         : graphic-user interface (GUI) version of the run-through script

DriverFalcon.m      : run-through script with provided examples, also possible to integrate more example and to       perform    manual tuning
