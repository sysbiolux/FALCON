% ======================================
% Driver script (run-through) for FALCON
% ======================================

% FalconInstall % In case the Falcon toolbox has not yet been added to Matlab's path
clc, clear all % clear screen and workspace 

% Choose your model example [1-4]
Model_Example = 4;

% 1 = Pipeline example
% 2 = PDGF model
% 3 = CellNOpt example
% 4 = Apoptosis model

% Define optmisation options
optRound=3; % Number of optimisation round
MaxFunEvals=3000; % Number of maximal function being evaluated (3000 = default)
MaxIter=3000; % Number of maximal iteration being evaluated (3000 = default)
Parallelisation=3; % Use multiple cores for optimisation? (0=no, 1=yes)
HLbound=0.5; % Qualitative threshold between high and low inputs
Forced=1; % Define whether single inputs and Boolean gates are forced to probability 1 
InitIC=2; % Initialise parameters' distribution (1=uniform, 2=normal)

% Define plotting and saving (0=no, 1=yes)
PlotFitEvolution    = 0; % Graph of optimise fitting cost over iteration
PlotFitSummary      = 0; % Graph of state values at steady-state versus measurements (all in 1)
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
KONodes_Analysis    = 1; % Node knock-out analysis

% ===================================================
% |||||||||||||||||||||||||||||||||||||||||||||||||||
% Click "Run" or press "F5" to start the optimisation
% |||||||||||||||||||||||||||||||||||||||||||||||||||
% ===================================================

%% Please modify the following part of the script for manual tuning

% Read model and measurement files 
if Model_Example == 1
    InputFile=[filesep 'ExampleDatasets' filesep 'example_model.txt'];
    MeasFile=[filesep 'ExampleDatasets' filesep 'example_meas.txt'];
elseif Model_Example == 2
    InputFile=[filesep 'ExampleDatasets' filesep 'PDGF_model.xlsx'];
    MeasFile=[filesep 'ExampleDatasets' filesep 'PDGF_meas.xlsx'];
elseif Model_Example == 3
    InputFile=[filesep 'ExampleDatasets' filesep 'CNO_model.xlsx'];
    MeasFile=[filesep 'ExampleDatasets' filesep 'CNO_data.xlsx'];
elseif Model_Example == 4
    InputFile=[filesep 'ExampleDatasets' filesep 'Apoptosis_model.xlsx'];
    MeasFile=[filesep 'ExampleDatasets' filesep 'Apoptosis_meas.xlsx'];
else
    error('Please specifiy the range of number from 1 to 4')
end

% Create a save folder
SaveFolderName=['Results_' datestr(now)];
FinalFolderName=strrep(SaveFolderName, ':', '.'); 
mkdir(FinalFolderName) % Automatically generate a folder for saving

% Build a FALCON model for optimisation
estim=FalconMakeModel([InputFile],[MeasFile],HLbound,Forced); %make the model

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

% Optimization
toc_all=[];
x_all=[];
fval_all=[];

if Parallelisation == 1

    parfor counter=1:optRound %'parfor' will run the parallel computing toolbox
        tic
        k=FalconIC(estim,IC_Dist); %initial conditions
        [xval,fval]=FalconObjFun(estim,k); %objective function
        Fitting_Time=toc;
        %collecting the run's outcome
        toc_all=[toc_all; Fitting_Time];
        x_all=[x_all; xval];
        fval_all=[fval_all; fval];
    end
 
else
    h = waitbar(0,'Please wait...');
    for counter=1:optRound %'parfor' will run the parallel computing toolbox
        waitbar(counter/optRound,h,sprintf('Running Optimisation Round %d out of %d ...',counter,optRound))
        tic
        k=FalconIC(estim,IC_Dist); %initial conditions
        [xval,fval]=FalconObjFun(estim,k); %objective function
        Fitting_Time=toc;
        %collecting the run's outcome
        toc_all=[toc_all; Fitting_Time];
        x_all=[x_all; xval];
        fval_all=[fval_all; fval];
    end
    close(h);
end

fxt_all=[fval_all x_all toc_all];

beep; pause(0.5); beep;

%Retrieving the results
[bestx,meanx,stdx]=FalconResults(fxt_all,estim.param_vector,FinalFolderName);

estim.MaxTime = mean(fxt_all(:,end))*3;
estim.Results.Optimisation.FittingCost = fxt_all(:,1);
estim.Results.Optimisation.FittingTime = fxt_all(:,end);
estim.Results.Optimisation.ParamNames = estim.param_vector;
estim.Results.Optimisation.BestParams = bestx;
estim.Results.Optimisation.StateNames = estim.state_names;

% Analyzing the evolution of fitting cost
if PlotFitEvolution == 1
  estim=FalconFitEvol(estim,IC_Dist,FinalFolderName);
end

%% Re-simulate results based on the best optimised parameter set
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
estim.MeanStateValueAll=MeanStateValueAll; estim.bestx=bestx;

%% Displayed optimised network with weights
if PlotBiograph == 1
    FalconShowNetwork(estim, PlotAllBiographs, FinalFolderName)
end

%% Obtaining robust estimates for the parameters by resampling
if Resampling_Analysis == 1;
    optRound_Sampling=optRound;
    CV_cutoff = 10; % Percent cut-off for large coefficient of variation
    [~,~,estim]=FalconResample(estim, bestx, optRound_Sampling, NDatasets, CV_cutoff, FinalFolderName);
end

%% Sensitivity analysis
if LPSA_Analysis == 1 
    optRound_LPSA=1;
    if Fast_Option == 1
        IsFast='fast';
    else
        IsFast='slow';
    end
    Estimated_Time_LPSA=mean(fxt_all(:,end))*optRound_LPSA*LPSA_Increments*3*length(estim.param_vector);
    disp(['Estimated Time for LPSA analysis (fast): ' num2str(Estimated_Time_LPSA) ' seconds']); beep; pause(3); beep; 
    [~, estim]=FalconLPSA(estim, bestx, MeasFile, HLbound, optRound_LPSA, LPSA_Increments, IsFast, Parallelisation, FinalFolderName);
end

%% Knock-out analysis
if KO_Analysis == 1;
    optRound_KO=1;
    Estimated_Time_KO=mean(fxt_all(:,end))*optRound_KO*length(estim.param_vector);
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep; 
    estim=FalconKO(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO, FinalFolderName);
end

%% Node Knock-out analysis
if KONodes_Analysis == 1;
    optRound_KO=1;
    Estimated_Time_KO_Nodes=mean(fxt_all(:,end))*optRound_KO*length(estim.NrStates);
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO_Nodes) ' seconds']); beep; pause(3); beep; 
    estim=FalconKONodes(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO,FinalFolderName);
end
%% Guided to display results in estim.Results
disp(' ')
disp('================================================')
disp('The optimisation and analyses are completed')
disp('================================================')

disp(' ')
disp('Please also check the results in "estim.Results"')
disp(estim.Results)
save([pwd filesep FinalFolderName filesep 'estim_Results'])
disp('================================================')

%% Model Prediction

if Model_Example == 2
    ValidationDataset = [filesep 'ExampleDatasets' filesep 'PDGF_validation.xlsx'];
    FalconValidation(estim,bestx,ValidationDataset,FinalFolderName);
end

%% Re-plot the figures (in case the plotting crashed during the analysis)
% FalconPlots(estim);

% === End of the script === %
