% ======================================
% Driver script (run-through) for FALCON
% ======================================
%%%
% FalconInstall % In case the Falcon toolbox has not yet been added to Matlab's path
clc, clear all % clear screen and workspace 

% Choose your model example [1-4]
Model_Example = 3;

% 1 = Pipeline example
% 2 = PDGF model
% 3 = CellNOpt example
% 4 = Apoptosis model

% Define optmisation options
optRound=10; % Number of optimisation round
MaxFunEvals=3000; % Number of maximal function being evaluated (3000 = default)
MaxIter=3000; % Number of maximal iteration being evaluated (3000 = default)
Parallelisation=1; % Use multiple cores for optimisation? (0=no, 1=yes)
HLbound=0.5; % Qualitative threshold between high and low inputs
InitIC=2; % Initialise parameters' distribution (1=uniform, 2=normal)

% Define plotting and saving (0=no, 1=yes)
PlotFitEvolution    = 0; % Graph of optimise fitting cost over iteration
PlotFitSummary      = 1; % Graph of state values at steady-state versus measurements (all in 1)
PlotFitIndividual   = 0; % Graph of state values at steady-state versus measurements (individual)
PlotHeatmapCost     = 0; % Heatmaps of optimal costs for each output for each condition absolute cost
PlotStateSummary    = 0; % Graph of only state values at steady-sate (all in 1)
PlotStateEvolution  = 0; % Graph of state values evolution over the course of the simulation (two graphs)
PlotBiograph        = 0; % Graph of network topology, nodes activities, and optimised parameters
PlotAllBiographs    = 0; % (Only for machines with strong GPUs) Plot all Biographs above

% Additional analyses after the optimisation with the default setting (0=no, 1=yes)
Resampling_Analysis = 0; % Resampling of experimental data and re-optimise
NDatasets           = 10;% Number of artificial datasets from which to resample.

LPSA_Analysis       = 0; % Local parameter sensitivity analysis
Fast_Option         = 1; % Performing faster LPSA by stopping if fitting costs go over a set threshold value
LPSA_Increments     = 4; % Number of increments for LPSA. Increase for finer resolution

KO_Analysis         = 0; % Parameter knock-out analysis

KO_Nodes_Analysis_eff = 1; % test different KO efficencies on each node and analyse the entire network based on this information
efficency_range = [1]; % indicate the different KO efficencies you want to test (vector from 0 to 1)
% ===================================================
% |||||||||||||||||||||||||||||||||||||||||||||||||||
% Click "Run" or press "F5" to start the optimisation
% |||||||||||||||||||||||||||||||||||||||||||||||||||
% ===================================================

%% Please modify the following part of the script for manual tuning

FALCONFolder = fileparts(which('DriverFalcon'));

% Read model and measurement files 
if Model_Example == 1
    InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'example' filesep 'example_model.txt'];
    MeasFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'example' filesep 'example_meas.txt'];
elseif Model_Example == 2
    InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'PDGF' filesep 'PDGF_model.xlsx'];
    MeasFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'PDGF' filesep 'PDGF_meas.xlsx'];
elseif Model_Example == 3
 InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'Greta' filesep 'crosstalk_model_final_IKBA.xlsx'];
    MeasFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'Greta' filesep 'crosstalk_parental_16h.xlsx'];
elseif Model_Example == 4
    InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'Greta' filesep 'crosstalk_model_final_IKBA.xlsx'];
    MeasFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'Greta' filesep 'crosstalk_izi_cond_16h.xlsx'];
else
    error('Please specifiy the range of number from 1 to 4')
end

% Create a save folder
SaveFolderName=['Results_' datestr(now)];
FinalFolderName=strrep(SaveFolderName, ':', '.'); 
mkdir(FinalFolderName) % Automatically generate a folder for saving

% Build a FALCON model for optimisation
estim=FalconMakeModel(InputFile,MeasFile,HLbound); %make the model

% Define optimisation options
estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
estim.SSthresh=eps;
if InitIC == 1
    IC_Dist ='uniform';
elseif InitIC == 2
    IC_Dist ='normal';
else
    error('Please choose the initial parameter distribution')
end




%% Re-simulate results based on the parental cell line

param_par_ikba= []
param_par_ikba= [0.013;0.500;0.379;0.354;0.981;0.071;0.736;0.006;0.921;0.679;0.004;0.731;0.022;0.006;0.019;0.978; ...
0.008;0.012;0.005;0.500;0.473;0.963;0.850;0.523;0.264;0.321;0.025;0.269;0.013];
 
%parental
param_par_ikba= [param_par_ikba; param_par_ikba(14)];
param_par_ikba(14) = [];
[MeanStateValueAll_Parental, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,param_par_ikba,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

%ikba superrepressor
param_par_ikba2= param_par_ikba;
param_par_ikba2(29)=1 ;
[MeanStateValueAll_Parental_ikba, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,param_par_ikba2,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

%XIAP kO   
param_par_xiap= param_par_ikba;
param_par_xiap(28) = 0;
[MeanStateValueAll_Parental_XIAP, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,param_par_xiap,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);





%% Re-simulate results based on the conditioned cell line

param_cond_ikba= []
param_cond_ikba= [0.013; 0.501; 0.004; 0.365;0.984;0.192;0.736;0.005;0.795;0.680;0.008;0.730;0.021;0.006;0.016; ...
            0.979;0.013;0.020;0.005;0.499;0.118;0.954;0.859;0.874;0.264;0.320;0.026;0.270;0.871];
%condtioned     
param_cond_ikba= [param_cond_ikba; param_cond_ikba(14)];
param_cond_ikba(14) = [];
[MeanStateValueAll_cond, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,param_cond_ikba,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
%conditioned super repressor
param_cond_xiap= param_cond_ikba;
param_cond_xiap(29)=1 ;
[MeanStateValueAll_cond_ikba, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,param_cond_xiap,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

% conditioned XIAP KO
param_par_xiap= param_par_ikba;
param_par_xiap(28) = 0;
[MeanStateValueAll_cond_XIAP, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,param_par_xiap,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
