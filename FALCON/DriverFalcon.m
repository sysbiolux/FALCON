% ======================================
% Driver script (run-through) for FALCON
% ======================================
%%%
% FalconInstall % In case the Falcon toolbox has not yet been added to Matlab's path
clc, clear all % clear screen and workspace 

% Choose your model example [1-4]
Model_Example = 2;

% 1 = Pipeline example
% 2 = PDGF model
% 3 = CellNOpt example
% 4 = Apoptosis model

% Define optmisation options
optRound=4; % Number of optimisation round
MaxFunEvals=3000; % Number of maximal function being evaluated (3000 = default)
MaxIter=3000; % Number of maximal iteration being evaluated (3000 = default)
Parallelisation=1; % Use multiple cores for optimisation? (0=no, 1=yes)
HLbound=0.5; % Qualitative threshold between high and low inputs
InitIC=2; % Initialise parameters' distribution (1=uniform, 2=normal)

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
NDatasets           = 10;% Number of artificial datasets from which to resample.

LPSA_Analysis       = 1; % Local parameter sensitivity analysis
Fast_Option         = 0; % Performing faster LPSA by stopping if fitting costs go over a set threshold value
LPSA_Increments     = 4; % Number of increments for LPSA. Increase for finer resolution

KO_Analysis         = 1; % Parameter knock-out analysis

KO_Nodes_Analysis_eff = 1; % test different KO efficencies on each node and analyse the entire network based on this information
efficency_range = [0:0.5:1]; % indicate the different KO efficencies you want to test (vector from 0 to 1)
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
    %%%for xls
    InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'PDGF' filesep 'PDGF_model.xlsx'];
    %%%for sif
%     InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'PDGF' filesep 'PDGF_model.sif'];
    %%%for xls
    MeasFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'PDGF' filesep 'PDGF_meas.xlsx'];
    %%%for csv
    MeasFile={'PDGF_meas_in.csv';'PDGF_meas_out.csv';'PDGF_meas_err.csv'};
elseif Model_Example == 3
    InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'CNO' filesep 'CNO_model.xlsx'];
    MeasFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'CNO' filesep 'CNO_data.xlsx'];
elseif Model_Example == 4
    InputFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'Apoptosis' filesep 'Apoptosis_model.xlsx'];
    MeasFile=[FALCONFolder filesep 'ExampleDatasets' filesep 'Apoptosis' filesep 'Apoptosis_meas.xlsx'];
else
    error('Please specifiy the range of number from 1 to 4')
end

% Create a save folder
SaveFolderName=['Results_' datestr(now)]; FinalFolderName=strrep(SaveFolderName, ':', '.'); 
mkdir(FinalFolderName) % Automatically generate a folder for saving

% Build a FALCON model for optimisation
estim=FalconMakeModel(InputFile,MeasFile,HLbound); %make the model

% Define optimisation options
estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
estim.SSthresh=eps;
if InitIC == 1,     IC_Dist ='uniform'; elseif InitIC == 2, IC_Dist ='normal'; else,    error('Please choose the initial parameter distribution'), end

% estim.Lambda=0.5;
% estim.Reg='L1';



% Optimization
toc_all=[]; x_all=[]; fval_all=[];

if Parallelisation
    parfor counter=1:optRound %'parfor' will run the parallel computing toolbox
        tic, k=FalconIC(estim,IC_Dist); %initial conditions
        
        [xval,fval]=FalconObjFun(estim,k); %objective function
        
        toc_all=[toc_all; toc]; x_all=[x_all; xval]; fval_all=[fval_all; fval];
    end

else
    for counter=1:optRound %'parfor' will run the parallel computing toolbox
        tic, k=FalconIC(estim,IC_Dist); %initial conditions
        
        [xval,fval]=FalconObjFun(estim,k); %objective function
        
        toc_all=[toc_all; toc]; x_all=[x_all; xval]; fval_all=[fval_all; fval];
    end
    close(h);
end

fxt_all=[fval_all x_all toc_all];

beep; pause(0.5); beep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Retrieving the results

[bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);

estim.MaxTime = mean(fxt_all(:,end))*3;
estim.Results.Optimisation.FittingCost = fxt_all(:,1);
estim.Results.Optimisation.FittingTime = fxt_all(:,end);
estim.Results.Optimisation.ParamNames = estim.param_vector;
estim.Results.Optimisation.BestParams = bestx;
estim.Results.Optimisation.StateNames = estim.state_names;

%%% Re-simulate results based on the best optimised parameter set
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
estim.MeanStateValueAll=MeanStateValueAll; estim.bestx=bestx;

save([FinalFolderName filesep 'OptimizedModel'])

%%
%%% Analyzing the evolution of fitting cost
if PlotFitEvolution
  estim=FalconFitEvol(estim,IC_Dist,FinalFolderName);
end

%%% Displayed optimised network with weights
if PlotBiograph
    FalconShowNetwork(estim, PlotAllBiographs, FinalFolderName)
end

%%% Obtaining robust estimates for the parameters by resampling
if Resampling_Analysis
    optRound_Sampling=optRound;
    CV_cutoff = 10; % Percent cut-off for large coefficient of variation
    [~,~,estim]=FalconResample(estim, bestx, optRound_Sampling, NDatasets, CV_cutoff, FinalFolderName);
end

%%% Sensitivity analysis
if LPSA_Analysis 
    optRound_LPSA=1;
    if Fast_Option
        IsFast='fast';
    else
        IsFast='slow';
    end
    Estimated_Time_LPSA=mean(fxt_all(:,end))*optRound_LPSA*LPSA_Increments*3*length(estim.param_vector);
    disp(['Estimated Time for LPSA analysis (fast): ' num2str(Estimated_Time_LPSA) ' seconds']); beep; pause(3); beep; 
    [~, estim]=FalconLPSA(estim, bestx, MeasFile, HLbound, optRound_LPSA, LPSA_Increments, IsFast, Parallelisation, FinalFolderName);
end

%%% Knock-out analysis
if KO_Analysis
    optRound_KO=1;
    Estimated_Time_KO=mean(fxt_all(:,end))*optRound_KO*length(estim.param_vector);
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep; 
    estim=FalconKO(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO, FinalFolderName);
end

%%% Nodes Knock-out analysis
% if KO_Nodes_Analysis
%     optRound_KO=1;
%     Estimated_Time_KO=mean(fxt_all(:,end))*optRound_KO*(length(estim.state_names)-length(estim.Input_idx(1,:)));
%     disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep; 
%     estim=FalconKONodes(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO, FinalFolderName);
% end

%%% Nodes Knock-out efficency analysis
if KO_Nodes_Analysis_eff
    optRound_KO=1;
    big_matrix = []
    
    for j= 1:length(efficency_range)
        estim.efficency= num2str(efficency_range(j))
        Estimated_Time_KO=mean(fxt_all(:,end))*optRound_KO*(length(estim.state_names)-length(estim.Input_idx(1,:)));
        disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep;
        estim=FalconKONodes_eff(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO, FinalFolderName);
        big_matrix (:,:,j)= estim.Results.KnockOutNodes.StateValue_screen;  
    end
   
    
    % % single figures per knock-out
    Nodes=estim.state_names;
    Nodes(estim.Input_idx(1,:))=[];
    for counter1 = 1: length(Nodes)
%         knock_out = estim.state_names(counter1)
%         matrix_split = squeeze(big_matrix(:,counter,:))
        var = ((counter1 - 1) * size(estim.Input,1)) + 1
         figure;
         h= suptitle(['Effect of ', char(Nodes(counter1)), ' Knock-out on:'])
         set(h,'FontSize',20,'FontWeight','bold')
    for counter = 1: length(estim.state_names)
        
        subplot(round(length(estim.state_names)/6), 6, counter)
        protein = estim.state_names(counter)
        matrix = squeeze(big_matrix(var : var + size(estim.Input,1)-1,counter,:))
        imagesc(matrix)
        title(estim.state_names(counter))
        xlabel('perturbation')
        set(gca, 'XTick', 1:1:11)
        set(gca, 'XTickLabel', efficency_range)       
    end
    end
    
     % % figure per condition
    Nodes=estim.state_names;
    Nodes(estim.Input_idx(1,:))=[];
    for counter1 = 1: size(estim.Input,1)
%         knock_out = estim.state_names(counter1)
%         matrix_split = squeeze(big_matrix(:,counter,:))
        var = ((counter1 - 1) * length(Nodes)) + 1
         figure;
         h= suptitle(['Effect of condition', num2str((counter1)), ' Knock-out on:'])
         set(h,'FontSize',20,'FontWeight','bold')
    for counter = 1: length(estim.state_names)
        subplot(round(length(estim.state_names)/6), 6, counter)
        protein = estim.state_names(counter)
        matrix = squeeze(big_matrix(counter1 : size(estim.Input,1) : end, counter,:))
        imagesc(matrix)
        title(estim.state_names(counter))
        xlabel('perturbation')
        ylabel ('Knock-out') 
        set(gca, 'XTick', 1:1:11)
        set(gca, 'XTickLabel', efficency_range)   
        set(gca, 'YTick', 1:length(Nodes))
        set(gca, 'YTickLabel', Nodes)  
    end
    end

end

%%% Guided to display results in estim.Results
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
    ValidationDataset = [filesep 'ExampleDatasets' filesep 'PDGF' filesep 'PDGF_validation.xlsx'];
    FalconValidation(estim,bestx,ValidationDataset,FinalFolderName);
end

%% Re-plot the figures (in case the plotting crashed during the analysis)
FalconPlots(estim);

% === End of the script === %
