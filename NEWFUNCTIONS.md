# New FALCON functions

This new version includes functions that are not documented in the original Bioinformatics paper. Hopefully this document provides enough information to guide the use of these functions in a meaningful way.

## Regularizations

Regularization is a technique to include additional constrains into an optimization problem. This is done by adding to the cost function, a function of the parameters themselves. For example, Tikhonov regularization, also called **L2** or **ridge regression** or **weight decay** considers the sum of the squares of the parameters values. In general, regularizations materialize prior knowledge about the parameter distributions, and so are very Bayesian in nature. For example, the Tikhonov penality induces smaller parameter weights, something that makes machine learning models usually more generalizable and moderates overfitting. By increasing incrementally the importance of this additional term in the objective function while controling for the total goodness of fit, it is possible to choose a threshold to determine the 'best' model. 
Some regularization schemes are now available for FALCON:

### L1/2

This introduces a partial norm on all parameters. This will induce a simplification in the topology of the network, in the sense that stronger regularization will force the optimization to *choose* between additive edges to the same nodes and will pick the most important of competitive interactions. It will also have a pruning effect on the negative edges.

### L1

This introduces a simple L1 norm on all parameters. Because of the Bayesian structure in FALCON (the sum of activating edges must be equal to 1) this is not useful for activating edges. This will however prune out less important inhibiting edges. 

### L2

This introduces a simple L2 norm on all parameters. Because of the Bayesian structure in FALCON (the sum of activating edges must be equal to 1) this is not useful for activating edges, actually pushing competitive edges to have more equal contributions. It is not clear what biological assumption this materializes. This will however decrease the weight of less important inhibiting edges. 

### L1Groups (global model)

This penalizes the sum of the absolute differences to the average parameter value for groups of parameters. It will make the parameter values for the same edge more equal across contexts, where contexts are for example cell lines or patients or treatments. By screening regularization strengths, this can determine which edges are probably context specific and which ones are probably not given the data. However in most cases we think this is not the best option when the number of contexts is more than 2. See LCluster below.

### L1Smooth (global model)

This penalizes the sum of cumulative differences across contexts, when the contexts have a natural order, for example timepoints, spatial location, increasing amount of genetic damage etc. This will 'smooth' parameter values and assume that differences between two successive contexts should be minimized. By screening regularization strengths, this can determine which edges are probably context specific (timepoint a change occurs for example) and which ones are probably not given the data. In spirit, this is similar to using convolutions for neural networks, and in general materializes the assumption of local signal correlation.

### LDrugs (global model) 

Combines L1/2 with LGroups. Not well tested.

### LTriple (global model)

Combines L1/2 with LGroups with LSmooth. Not well tested.

### LCluster (global model)

This is an original contribution of ours (see De Landtsheer _et al._, 2018, Frontiers on Physiology, xxx). By summing the local deviations to the expected density for all n-dimensional rectangles in the parameter space, we obtain a penality which we demonstrated is useful to group contexts together when there is a high number of them, and that therefore a certain level of heterogeneity is expected. By screening regularization strengths, this can determine which edges are probably context specific and which ones are probably not given the data.

### PruneCluster (global model)

Combines L1/2 with LCluster. Not well tested.

# How to use regularization with FALCON

See the driver script **Driver_Regularization.m**. It is not possible at the moment to use regularized optimizations with the GUI.
    - The variable *estim.Reg* should contain a string and defines the type of regularization applied. 
    - The variable *estim.RegMatrix* should contain for global models (when there are several contexts) a N-by-P matrix with parameter ID number for N edges and P contexts. For example, a 5-edge network across 3 contexts would basically have the following matrix: [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15]. It is automatically generated but can be customized before optimization, for example to exclude some edges from regularization or group parameters even further. In the above example we could want to have [1 2 3 4 5 6; 7 8 9 10 11 12] and not care about the last edge.
    - the variable *estim.Lambda* should contain a vector with the values of the lambda parameter, or strength of regularization. We recommend screening large ranges at first and increase sampling as much as computationally possible, for example choose values between 2^-20 to 1, with half-log increments (40 simulations).
    - For mixed regularizations, we need two or three vectors of lambda values (depending on the case), as we will perform a grid search in two or three dimensions.
    
After the optimizations are over, the one thing to look for is the most appropriate balance between goodness-of-fit (Mean Squared Error) and model size (number of non-null parameters). Different metrics can be used, we usually prefer the BIC, however AIC often leads to the same conclusion about model topology. We consider 0.01 as the threshold for a parameter to be different from zero.


## Data preprocessing

The function *FalconPreProcess.m* allows the manipulation of the measurement data in a number of ways. It needs to be run when the optimization problem has been created, before the optimization process.

### Normalization

*[estim2]=FalconPreProcess(estim, 'normalize', [a b])* will min-max normalize the values in the dataset between a and b, with a < b (typically 0-1). Normalization is performed linearly, independently for each analyte (each column in the dataset).

### Bootstrapping

*[estim2]=FalconPreProcess(estim, 'bootstrap', S)* where S is a string with possible values 'rows', 'columns', or 'both' will perform bootstrap sampling (resample with replacement to the original dataset size). This is a simple way to assess the confidence in parameter values estimates given the amount of signal in the data, by looking at the parameter value distributions over many bootstraps.

### Randomization

*[estim2]=FalconPreProcess(estim, 'randomize', p)* where p is a scalar with 0 < p < 1 will perform random permutations on a fraction p of the datapoints. This is an easy way to assess the influence of the data on the model optimization. Reaching the same conclusion with the original and randomized dataset is a sign that the conclusions are not well supported.

### Subsampling

*[estim2]=FalconPreProcess(estim, 'subsample', p)* where p is a positive scalar will subsample the original dataset to a fraction p of its original size. If p < 1 the resulting dataset will be smaller than the original one (undersampling, useful to speed up computations on large datasets). If p > 1 the result is oversampling, which is a less aggressive way to perturbate the dataset than bootstrapping, as all original datapoints are still included, and therefore better in the case of very small datasets.

### Noise addition

*[estim2]=FalconPreProcess(estim, 'noise', p)* where p is a positive scalar will add a fraction p of the error on the measurements to the measurement values, de facto generating one of the possible instances of measurement values given their distribution. This is yet another way to perturb the data, in a way that preserves the error matrix. This is useful when some measurement have a much higher standard error than others.

It is possible to combine different layers of preprocessing in one function call, like in: *[estim2]=FalconPreProcess(estim, 'normalize', [-0.2 0.9], 'bootstrap', 'both','randomize', 0.5, 'subsample', 0.5)*

