%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      _______  _______  _        _______  _______  _          %
%     (  ____ \(  ___  )( \      (  ____ \(  ___  )( (    /|   %
%     | (    \/| (   ) || (      | (    \/| (   ) ||  \  ( |   %
%     | (__    | (___) || |      | |      | |   | ||   \ | |   %
%     |  __)   |  ___  || |      | |      | |   | || (\ \) |   %
%     | (      | (   ) || |      | |      | |   | || | \   |   %
%     | )      | )   ( || (____/\| (____/\| (___) || )  \  |   %
%     |/       |/     \|(_______/(_______/(_______)|/    )_)   %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Current version: v1.2                      %
%                      (February 2019)                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                              
%    CONTEXTUALIZATION OF DBN MODELS OF BIOLOGICAL NETWORKS    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Sébastien De Landtsheer: sebastien.delandtsheer@uni.lu      %
%                           seb@delandtsheer.com               %
%  Prof Thomas Sauter:      thomas.sauter@uni.lu               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          University of Luxembourg 2019 GPLv3                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cite: De Landtsheer et al., 2017, Bioinformatics, 33, 3431–6 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all % clear screen and workspace 


% ===========================
% Define optimization options
% ===========================
optRound = 4; % Number of optimisation round
MaxFunEvals = 3000; % Number of maximal function being evaluated (3000 = default)
MaxIter = 3000; % Number of maximal iteration being evaluated (3000 = default)
Parallelisation = 1; % Use multiple cores for optimisation? (0=no, 1=yes)
HLbound = 0.5; % Qualitative threshold between high and low inputs
InitIC = 2; % Initialise parameters' distribution (1=uniform, 2=normal, 3=zeros)

% Define plotting and saving (0=no, 1=yes)
PlotFitEvolution    = 1; % Graph of optimise fitting cost over iteration
PlotFitSummary      = 1; % Graph of state values at steady-state versus measurements (all in 1)
PlotFitIndividual   = 0; % Graph of state values at steady-state versus measurements (individual)
PlotHeatmapCost     = 1; % Heatmaps of optimal costs for each output for each condition absolute cost
PlotStateSummary    = 1; % Graph of only state values at steady-sate (all in 1)
PlotStateEvolution  = 1; % Graph of state values evolution over the course of the simulation (two graphs)
PlotBiograph        = 0; % Graph of network topology, nodes activities, and optimised parameters
PlotAllBiographs    = 0; % (Only for machines with strong GPUs) Plot all Biographs above

% Additional analyses after the optimisation with the default setting (0=no, 1=yes)
Resampling_Analysis = 1; % Resampling of experimental data and re-optimise
NDatasets           = 50;% Number of artificial datasets from which to resample.

LPSA_Analysis       = 1; % Local parameter sensitivity analysis
Fast_Option         = 0; % Performing faster LPSA by stopping if fitting costs go over a set threshold value
optRound_LPSA       = 5; % Number of optimizations for each perturbed datapoint
LPSA_Increments     = 3; % Number of increments for LPSA. Increase for finer resolution

KO_Analysis         = 1; % Parameter knock-out analysis
KO_Nodes_Analysis   = 1; % Node knock-out analysis
optRound_KO         = 5; % Number of optimizations for each KO datapoint. 
KO_Nodes_fast       = 1; % 

KO_Nodes_Analysis_eff = 1; % test different KO efficencies on each node and analyse the entire network based on this information
efficiency_range = 0:0.1:1; % indicate the different KO efficencies you want to test (vector from 0 to 1)

% Model and measurement files
InputFile = 'PDGF_model.xlsx'; %for xls
% InputFile='PDGF_model.sif'; %for sif
MeasFile = 'PDGF_meas.xlsx'; %for xls
% MeasFile={'PDGF_meas_in.csv';'PDGF_meas_out.csv';'PDGF_meas_err.csv'}; %for csv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Creates a save folder
SaveFolderName = ['Results_' datestr(now)]; FinalFolderName=strrep(SaveFolderName, ':', '.'); 
mkdir(FinalFolderName) % Automatically generate a folder for saving

% Builds a FALCON model for optimisation
estim = FalconMakeModel(InputFile, MeasFile, HLbound); %make the model

% Defines optimisation options
estim.options = optimoptions('fmincon', 'TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-10, 'MaxFunEvals', MaxFunEvals, 'MaxIter', MaxIter); % Default setting
estim.SSthresh=eps;
if InitIC == 1, IC_Dist = 'uniform'; elseif InitIC == 2, IC_Dist = 'normal'; elseif InitIC == 3, IC_Dist = 'scratch'; else,    error('Please choose the initial parameter distribution'), end

% Optimization
toc_all = []; x_all = []; fval_all = [];

if Parallelisation
    parfor counter = 1:optRound %'parfor' will run the parallel computing toolbox
        tic, k = FalconIC(estim, IC_Dist); %initial conditions
        [xval, fval] = FalconObjFun(estim, k); %objective function
        toc_all = [toc_all; toc]; x_all = [x_all; xval]; fval_all = [fval_all; fval];
    end

else
    for counter = 1:optRound %'parfor' will run the parallel computing toolbox
        tic, k = FalconIC(estim, IC_Dist); %initial conditions
        [xval, fval] = FalconObjFun(estim, k); %objective function
        toc_all = [toc_all; toc]; x_all = [x_all; xval]; fval_all = [fval_all; fval];
    end
    close(h);
end

fxt_all = [fval_all x_all toc_all];
beep; pause(0.5); beep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Retrieving the results

[bestx, meanx, stdx, estim] = FalconResults(estim, fxt_all, estim.param_vector, FinalFolderName);

%%% Re-simulate results based on the best optimised parameter set
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim, bestx, [PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution], FinalFolderName);

save([FinalFolderName filesep 'OptimizedModel'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analyzing the evolution of fitting cost
if PlotFitEvolution
  estim = FalconFitEvol(estim, IC_Dist, FinalFolderName);
end

%%% Displayed optimised network with weights
if PlotBiograph
    FalconShowNetwork(estim, PlotAllBiographs, FinalFolderName)
end

%%% Obtaining robust estimates for the parameters by resampling
if Resampling_Analysis
    optRound_Sampling = optRound;
    CV_cutoff = 10; % Percent cut-off for large coefficient of variation
    [~,~,estim] = FalconResample(estim, bestx, NDatasets, 1, CV_cutoff, FinalFolderName);
end

%%% Sensitivity analysis
if LPSA_Analysis 
    if Fast_Option, IsFast = 'fast'; else, IsFast = 'slow'; end
    Estimated_Time_LPSA = mean(fxt_all(:, end)) * optRound_LPSA * LPSA_Increments * 3 * length(estim.param_vector);
    disp(['Estimated Time for LPSA analysis (fast): ' num2str(Estimated_Time_LPSA) ' seconds']); beep; pause(3); beep; 
    [~, estim] = FalconLPSA(estim, bestx, MeasFile, HLbound, optRound_LPSA, LPSA_Increments, IsFast, Parallelisation, FinalFolderName);
end

%%% Knock-out analysis
if KO_Analysis
    Estimated_Time_KO = mean(fxt_all(:, end)) * optRound_KO * length(estim.param_vector);
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep; 
    estim = FalconKO(estim, fxt_all, HLbound, optRound_KO, Parallelisation, FinalFolderName);
end

%%% Nodes Knock-out analysis
if KO_Nodes_Analysis
    Estimated_Time_KO = mean(fxt_all(:, end)) * optRound_KO * (length(estim.state_names) - length(estim.Input_idx(1,:)));
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep; 
    estim = FalconKONodes(estim, fxt_all, HLbound, optRound_KO, Parallelisation, FinalFolderName);
end

%%% Nodes Knock-out efficency analysis
if KO_Nodes_Analysis_eff
    Estimated_Time_KO = mean(fxt_all(:,end)) * numel(efficiency_range) * optRound_KO * (length(estim.state_names)-length(estim.Input_idx(1,:)));
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep;
    estim=FalconKONodes_eff(estim, fxt_all, efficiency_range, HLbound, optRound_KO, FinalFolderName);    

end

%%% Fast KO
if KO_Nodes_fast
    estim = FalconKONodes_fast(estim, fxt_all, HLbound, optRound_KO, Parallelisation, FinalFolderName);  
end

%%% Discriminate fast effects from long-term effects
X = estim.Results.KnockOutNodes.BIC_values;
Y = estim.Results.KnockOutNodesFast.BIC_values;
M1 = min(X,Y); M2 = min(X,Y);
figure, plot(X,Y,'.k', 'MarkerSize', 15), hold on
plot([M1, M2], [M1, M2], '-k');


% === End of the script === %
