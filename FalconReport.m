%% FALCON Network optimization report
%%
% Report created by Philippe at:
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
%%
% 
% <<Biograph_Exp1.jpg>>
% 
%% Network optimzation results
% Time the optimzation, best fitting cost and the weights of the
% optimized parameters. For each result the best as well as the mean and
% the standard deviation (S.D.) of the results is given.
% 
FalconResults(fxt_all,estim.param_vector,FinalFolderName)
%%
% The graph displays the goodness of the optimization for each protein.
% Green dots and bars represent the experimental data and the blue stars
% correspond to the optimized values. The numbers on the x-axis correspond
% to the experimental condition.
% 
% <<Fitting_plot.jpg>>
% 
%% 
% The heatmap displays all the optimized conditions gainst all the measured
% species. The color-gradient of the heatmap indicates the goodness of the fit for the different conditions.
% Conditions with a lighter color correspond to worse fitting costs.
% 
% <<CrossErrorHeatMap.jpg>>
% 

%% Resampling results
% FalconResample creates a new dataset based on the mean and standard deviation of the original dataset. 
% Useful to assess the variability of the parameter estimates, in
% conjunction with FalconLPSA.
% The resampled parameter values basd on the mean and the standard
% deviation are the following:
if isfield(estim.Results,'Resampling')
Parameter_names = estim.Results.Resampling.Parameters;
Resampling_results = estim.Results.Resampling.OptimisedParameter;
Resampling_SD = estim.Results.Resampling.OptimisedSD;

Heading=cell(1,3);
Heading(1,1)={'parameters'};
Heading(1,2)={'mean'};
Heading(1,3)={'S.D.'};

Resampling=[Parameter_names' num2cell(mean(Resampling_results)') num2cell(Resampling_SD')];

disp('Summary of resampling:')
disp(' ')
disp([Heading;Resampling])
disp(' ')
end
%%
% The figures correspond to the new artificial datasets created based on
% the experimental data.
% 
% <<Resampling.jpg>>
% 

%% Local parameter sensitivity results
% FalconLPSA performs a local parameters sensitivity analysis (LPSA) on an optimized model
if isfield(estim.Results,'LPSA')
   Parameter_names = estim.Results.LPSA.ParamNames;
   Resampling_cutoff = estim.Results.LPSA.CutOff
   LPSA_interpretation = estim.Results.LPSA.Interpretation
   LPSA_identifability = estim.Results.LPSA.Identifiability;
   
   
Heading=cell(1,2);
Heading(1,1)={'parameters'};
Heading(1,2)={'Identifiable?'};


LPSA=[Parameter_names' num2cell(LPSA_identifability') ];

disp('Summary of local parameter sensistivity analysis:')
disp(' ')
disp([Heading;LPSA])
disp(' ')
   
   
end
%%
% The figure displays the different parameters of the model. From the
% optimal parameter set (blue star), each parameter is perturbed in the
% direction of the upper and lower bound. The perturbed parameter value is
% fixed and all the other parameters are re-optimized. The red line
% corresponds to the cut-off calculated by the resampling.
% 
% <<Identifiability.jpg>>
% 
%% Knock-out results
% FalconKO creates new models for each knock-outed parameter by setting
% parameter value to 0 and re-optimise.
% The knock-out analysis identifies the parameters which are essential to
% correctly dexribe the experimental data. For this approach each model - 1
% parameter is compared against the reference model by the Akaike
% information criterion (AIC). If the removal of the parameter is not
% influencing the goodness of the fit, the AIC of the new model is the same
% or lower than the one in the refernce model. If the parameter is
% essential for the model structure, the AIC of the new model is increased
% compared t the reference model.
% 
if isfield(estim.Results,'KnockOut')
    Parameter_names = estim.Results.KnockOut.Parameters;
    AIC_values = estim.Results.KnockOut.AIC_values;
    KO_effect = estim.Results.KnockOut.KO_effect;
    KO_interpretation = estim.Results.KnockOut.Interpretation

    
Heading=cell(1,3);
Heading(1,1)={'parameters'};
Heading(1,2)={'AIC?'};
Heading(1,3)={'KO_effect'};

KO=[Parameter_names' num2cell(AIC_values') num2cell(KO_effect')];

disp('Summary of Knock out analysis:')
disp(' ')
disp([Heading;KO])
disp(' ')
    
end

%%
% The figure displays the obtained results for the knock-out analysis. The
% AIC for the reference model is colored in blue. Parameters having an
% effect on the model optimization after removal are displayed in black. 
% Parameters having no effect on the model optimization after removal are 
% displayed in green.
% 
% <<Knock-Outs.jpg>>
% 