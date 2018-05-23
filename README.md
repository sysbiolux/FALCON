**Contact info**: thomas.sauter@uni.lu or sebastien.delandtsheer@uni.lu  
**Please cite** : De Landtsheer _et al._, FALCON: A Toolbox for the Fast Contextualisation of Logical Networks (2017), Bioinformatics [Link](https://academic.oup.com/bioinformatics/article/33/21/3431/3897376)  
De Landtsheer _et al._, Using Regularization to Infer Cell Line Specificity in Logical Network Models of Signaling Pathways (2018) [Link](https://www.frontiersin.org/articles/10.3389/fphys.2018.00550/full)  
**Follow us on Twitter** : [FALCON toolbox](https://twitter.com/FALCON_toolbox)  
**Copyright**: This software is freely available for non-commercial users under the GPLv3 license.

# Version Info:

This updated version, released in October 2017, introduces new algorithms for systems analysis of regulatory networks and provides several bug fixes and improvements. Most notably, we removed several variables from the main functions that became irrelevant, got rid of nearly-mandatory Excel, made the plots more readable and intuitive, and added various regularization schemes to our objective function. We really hope these changes will improve the usability and usefulness of the toolbox and we encourage users to contribute to the project, and/or post their positive or negative experiences on the repository page.

Here is a list of the main changes:
* The 'Forced' variable has been removed and the related quantity fixed to 1. This makes the framework completely coherent from the mathematical point of view, as all nodes now are 'forced' to abide the rules of probability.
* The option to set 'high' and 'low' parameters, i.e. to bias the constrains in favor of some interactions based on prior knowledge becomes optional. We hope this makes the model input format clearer.
* Network model files and data files in the '.csv' format are now accepted. We hope this solves in an effective way issues related with the use of Excel on Mac, Linux, or older versions of Matlab. Example input files have been updated accordingly. As csv is the _de facto_ format for data science flat files, we might deprecate the '.txt' and Excel formats in later versions.
* Network models can be loaded in the '.sif' format. This makes FALCON model files readable by other software, for example CytoScape. Parameter names will be infered from the nodes names. All interactions will be considered linear (no Boolean gates).
* Data files now need to contain a first column termed 'annotation', which is a string (for example 'EGF+NoInhibitor-24h') associated with each condition, used for identifying this particular condition in later analyses and plots. Example input files have been updated accordingly.
* Each instance of FALCON will now generate a log file with basic info about data, model, hyperparameters and runtime info.
* The optimized metric is now MSE (Mean Squared Error). It is equal to the SSE (Sum of Squared Errors, used in the previous version) divided by the number of datapoints, and has the advantage of being a more sensible (not fool-proof though) measure to compare model fitness between different models and datasets.
* We spotted a bug in the use of Mersenne-Twister pseudo-random number generator, which seed is now really randomized according to the system status.
* The function _FalconKONodes_ has been expanded to include the option of partial knock-downs, in order to assess the effect of these on every other network node in the different conditions. It is meant to predict the effect of siRNA experiments, which rarely achieve 100% silencing.
* We introduce the new function _FalconPreProcess_ to handle data prior to optimization. Options include undersampling, oversampling, noise addition, normalization, bootstrapping and randomization along samples and/or conditions.
* The objective function _FalconObjFun_ now returns the MSE, SSE, AIC, BIC, and the number of parameters in the model (which could be different than the starting number in case of regularization).
* Many regularization schemes have been added: Straight: applying a partial norm to prune a network; Grouped: applying L1 norm on intra-group parameter variability across contexts to detect context-specific regulation; Smooth: applying L1 norm accross timepoints to detect timepoints where network re-wiring occurs; Uniformity: using local density disparity in the parameter space to cluster parameter together (see the related paper De Landtsheer _et. al._ (2018), Frontiers in Physiology, xxx). We also allow some combinations of several regularizations on the same problem. See the related documentation (Regularizations.md) for the use of regularized FALCON objective functions.
* In general, the toolbox has been brought to the field standard in term of annotation and documentation of open-source software.



***

# FALCON

> FALCON is a toolbox for the efficient contextualization of logical network models. Contextualized network models provide important qualitative and quantitative information about the system being modelled. Specifically, FALCON is well suited to assess the relative contributions of different signal transduction mechanisms to the behavior of the system at steady-state. The logical formulation of the interaction between different molecules is intuitive, and often is the only available information to systems biologists. A Bayesian interpretation of the logical ‘gates’ allows for an algebraic formulation of the system and an efficient calculation of the long-term steady-state of the system given the specified inputs. A gradient-descent optimizer (fmincon) is used to minimize the error function calculated as the sum of squared residuals between the simulated and experimentally measured nodes of interest. FALCON offers multiple types of analysis including identifiability analysis, systematic knock-out and differential regulation which can be applied to various biological applications.

***

# Example Applications

* [Toy model](https://github.com/sysbiolux/FALCON/wiki/Toy-example)
* [PDGF](https://github.com/sysbiolux/FALCON/wiki/PDGF) (Trairatphisan,P. _et al._ (2016))
* [MAPK](https://github.com/sysbiolux/FALCON/wiki/MAPK) (Saez-Rodriguez,J. _et al._ (2009))
* [Apoptosis](https://github.com/sysbiolux/FALCON/wiki/Apoptosis) (Schlatter,R. _et al._ (2009))

***

# Installation and requirements

In order to use the FALCON toolbox, the necessary scripts need to be included in the Matlab path. This can be done automatically by typing “FalconInstall” on the Matlab command window. The pipeline was developed under Matlab R2014b and has been successfully tested under Matlab R2015a and R2016a. It requires the build-in Matlab solver “fmincon”, which is a part of the Optimization Toolbox (http://nl.mathworks.com/products/optimization/), during the optimisation process. The parallel computing option in the FALCON pipeline which uses the “parfor” function requires the Parallel Computing Toolbox (http://nl.mathworks.com/help/distcomp/) and the plotting of network structure with optimised parameters requires the “biograph” function integrated in the Bioinformatics toolbox (http://nl.mathworks.com/help/bioinfo/) to be present.  

In case you encounter any problems while running the FALCON pipeline on earlier of later versions of Matlab, please do not hesitate to contact the authors or post an issue under: https://github.com/sysbiolux/FALCON/issues.

***

# Running the Falcon pipeline
* [Input files](https://github.com/sysbiolux/FALCON/wiki/Input-files)
* [Falcon functions](https://github.com/sysbiolux/FALCON/wiki/The-FALCON-functions)
* [Running Falcon](https://github.com/sysbiolux/FALCON/wiki/Running-FALCON)

***

# All variables
[All variables explained](https://github.com/sysbiolux/FALCON/wiki/Variables)

***

# Troubleshooting
The FALCON pipeline was rigorously tested before its release into the modelling community. Nevertheless, users might still encounter a few issues in the FALCON pipeline and we provide some potential solutions in this section.

* [The FALCON scripts and examples are not added to the path](https://github.com/sysbiolux/FALCON/wiki/The-FALCON-scripts-and-examples-are-not-added-to-the-path)
* [The FALCON pipeline does not accept certain logical gate combinations](https://github.com/sysbiolux/FALCON/wiki/The-FALCON-pipeline-does-not-accept-certain-logical-gate-combinations)
* [The FALCON pipeline crashed during the analytical process](https://github.com/sysbiolux/FALCON/wiki/The-FALCON-pipeline-crashed-during-the-analytical-process)

***

# Bibliography
* Lommel,M.J. _et al._ (2016) L-plastin Ser5 phosphorylation in breast cancer cells and in vitro is mediated by RSK downstream of the ERK/MAPK pathway. FASEB J., 30, 1218–33. 
[Link](http://www.fasebj.org/content/30/3/1218.long) 
* Saez-Rodriguez,J. _et al._ (2009) Discrete logic modelling as a means to link protein signalling networks with functional analysis of mammalian signal transduction. Mol. Syst. Biol., 5, 331. 
[Link](http://msb.embopress.org/content/5/1/331.long)
* Schlatter,R. _et al._ (2009) ON/OFF and beyond--a boolean model of apoptosis. PLoS Comput. Biol., 5, e1000595. 
[Link](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000595   )
* Trairatphisan,P. _et al._ (2016) A Probabilistic Boolean Network Approach for the Analysis of Cancer-Specific Signalling: A Case Study of Deregulated PDGF Signalling in GIST. PLoS One, 11, e0156223. 
[Link](http://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0156223  )
* Trairatphisan,P. _et al._ (2014) optPBN: an optimisation toolbox for probabilistic Boolean networks. PLoS One, 9, e98001.
[Link]( http://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0098001)
 


