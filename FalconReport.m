%% FALCON Network optimization report
%%
% Report performed by Philippe at:
datestr(clock,0)

%% Network structure
% FalconMakeModel creates the integrated optimization problem for Falcon from a list of interactions (InputFile)
% and a list of measured nodes in different conditions (MeasFile).
% The weights of interactions for single-input node(s) and constrains for the optimization problem are created.
% The following network structure was used for the optimization and
% contains the following nodes:
Network = estim.Interactions
Nr_of_nodes = estim.NrStates
Node_names = estim.state_names

%% Network optimzation results
% Time the optimzation, best fitting cost and the weights of the
% optimized parameters. For each result the best as well as the mean and
%  
% the standard deviation (S.D.) of the results is given.
FalconResults(fxt_all,estim.param_vector,FinalFolderName)
%%
% 
% <<Fitting_plot.jpg>>
% 
%%
% 
% <<CrossErrorHeatMap.jpg>>
% 

%% Resampling results
% FalconResample creates a new dataset based on the mean and standard deviation of the original dataset. 
% Useful to assess the variability of the parameter estimates, in
% conjunction with FalconLPSA.
% The resampled parameter values basd on the mean and the standard
% deviation are the following:
Parameter_names = estim.Results.Resampling.Parameters
Resampling_results = estim.Results.Resampling.OptimisedParameter
Resampling_SD = estim.Results.Resampling.OptimisedSD

%%
% 
% <<Resampling.jpg>>
% 
%% Local parameter sensitivity results
% FalconLPSA performs a local parameters sensitivity analysis (LPSA) on an optimized model
%%
% 
% <<Identifiability.jpg>>
% 
%% Knock-out results
% FalconKO creates new models for each knock-outed parameter by setting parameter value to 0 and re-optimise
%%
% 
% <<Knock-out.jpg>>
% 