FALCON
Fast algorithm for the contextualisation of logical networks.

Sébastien De Landtsheer, sebastien.delandtsheer@uni.lu

::: Contents :::

FalconInstall.m     : add Falcon scripts and example to the path (type FalconInstall on the Command Window)

DriverFalcon.m      : run-through script with provided examples, also possible to integrate more example and to perform manual tuning

FalconMakeModel.m   : creates the optimization problem from a list of interactions (InputFile) and a list of measured nodes in different conditions (MeasFile). The files can be .txt (see 'example_mod.txt' and 'example_meas.txt') or Excel (see 'PDGF.xlsx' and 'PDGF_meas.xlsx'). Automatically assigns the weights of interactions for single-input nodes and creates the constrains for the optimization problem.

FalconExpand        : automatic expansion of logical networks in order to facilitate rapid computation

FalconInt2File.m    : write the list of new interactions into a temporary file for re-importing (essential to generate perturbed model variants)

FalconIC.m          : automatically assigns initial parameter values which are drawn from a uniform or normal distribution.

FalconObjFun.m      : optimization with built-in objective function that simulates networks and returns the sum of squared errors between the simulated and expected values. Uses fmincon with the default algorithm (interior-point). Can be used inside a 'parfor' loop to perform multiple optimizations in parallel.

FalconObjFunPlus.m  : same as above but returns a global variable containing the costs over the course of the optimization. Useful for convergence diagnostics. Cannot be used inside a parfor loop.

FalconResults.m     : summarises the optimization output and display the best and average values for the parameters.

FalconFitEvol.m     : re-optimise the model and plot the evolution of fitting costs (3 rounds)

FalconSimul.m       : resimulates the optimized network and displays the summarised results as plots

FalconShowNetwork.m : uses the 'biograph' environment to display the interaction graphs of the networks.

FalconLPSA.m        : performs a local parameters sensitivity analysis on an optimized model.

FalconKO.m          : creates new models for each knock-outed parameter by setting parameter value to 0 and re-optimise

FalconResample.m    : creates a new dataset from the original one, given the mean and standard deviation of the original one. Useful to assess the variability of the parameter estimates, in conjonction with FalconLPSA.m

FalconValidation.m  : resimulates the optimized network and compare to the validation dataset

FalconPlots.m       : re-plot the figures from the chosen analysis (in case the plot(s) crashed during the analytical process)

FalconGUI.m         : graphic-user interface (GUI) version of the run-through script