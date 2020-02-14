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
%                      (July 2019)                         %
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

% Model and measurement files
InputFile = 'PDGF_model.xlsx'; %for xls
% InputFile='PDGF_model.sif'; %for sif
MeasFile = 'PDGF_meas.xlsx'; %for xls
% MeasFile={'PDGF_meas_in.csv';'PDGF_meas_out.csv';'PDGF_meas_err.csv'}; %for csv

% Builds a FALCON model for optimisation
estim = FalconMakeModel(InputFile, MeasFile); %make the model

%%%%
% modify optimization parameters here
%%%%

% Optimization
toc_all = []; x_all = []; fval_all = [];

if estim.Parallelisation
    parfor counter = 1:estim.optRound %'parfor' will run the parallel computing toolbox
        tic, k = FalconIC(estim); %initial conditions
        [xval, fval] = FalconObjFun(estim, k); %objective function
        toc_all = [toc_all; toc]; x_all = [x_all; xval]; fval_all = [fval_all; fval];
    end

else
    for counter = 1:estim.optRound %'parfor' will run the parallel computing toolbox
        tic, k = FalconIC(estim); %initial conditions
        [xval, fval] = FalconObjFun(estim, k); %objective function
        toc_all = [toc_all; toc]; x_all = [x_all; xval]; fval_all = [fval_all; fval];
    end
    close(h);
end

fxt_all = [fval_all x_all toc_all];
beep; pause(0.5); beep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Retrieving the results
[bestx, meanx, stdx, estim] = FalconResults(estim, fxt_all);

%%% Re-simulate results based on the best optimised parameter set
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim);

save([FinalFolderName filesep 'OptimizedModel'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analyzing the evolution of fitting cost
if estim.PlotFitEvolution
  estim = FalconFitEvol(estim);
end

%%% Displayed optimised network with weights
if estim.PlotBiograph
    FalconShowNetwork(estim)
end

%%% Obtaining robust estimates for the parameters by resampling
if estim.Resampling_Analysis
    optRound_Sampling = optRound;
    CV_cutoff = 10; % Percent cut-off for large coefficient of variation
    [~,~,estim] = FalconResample(estim, bestx, NDatasets, 3, CV_cutoff, FinalFolderName);
end

%%% Sensitivity analysis
if estim.LPSA_Analysis 
    if Fast_Option, IsFast = 'fast'; else, IsFast = 'slow'; end
    Estimated_Time_LPSA = mean(fxt_all(:, end)) * optRound_LPSA * LPSA_Increments * 3 * length(estim.param_vector);
    disp(['Estimated Time for LPSA analysis (fast): ' num2str(Estimated_Time_LPSA) ' seconds']); beep; pause(3); beep; 
    [~, estim] = FalconLPSA(estim, bestx, MeasFile, HLbound, optRound_LPSA, LPSA_Increments, IsFast, Parallelisation, FinalFolderName);
end

%%% Knock-out analysis
if estim.KO_Analysis
    Estimated_Time_KO = mean(fxt_all(:, end)) * optRound_KO * length(estim.param_vector);
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep; 
    estim = FalconKO(estim, fxt_all, HLbound, optRound_KO, Parallelisation, FinalFolderName);
end

%%% Nodes Knock-out analysis
if estim.KO_Nodes_Analysis
    Estimated_Time_KO = mean(fxt_all(:, end)) * optRound_KO * (length(estim.state_names) - length(estim.Input_idx(1,:)));
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep; 
    estim = FalconKONodes(estim, fxt_all, HLbound, optRound_KO, Parallelisation, FinalFolderName);
end

%%% Nodes Knock-out efficency analysis
if estim.KO_Nodes_Analysis_eff
    Estimated_Time_KO = mean(fxt_all(:,end)) * numel(efficiency_range) * optRound_KO * (length(estim.state_names)-length(estim.Input_idx(1,:)));
    disp(['Estimated Time for KO analysis: ' num2str(Estimated_Time_KO) ' seconds']); beep; pause(3); beep;
    estim=FalconKONodes_eff(estim, fxt_all, efficiency_range, HLbound, optRound_KO, FinalFolderName);    

end

%%% Fast KO
if estim.KO_Nodes_fast
    estim = FalconKONodes_fast(estim, fxt_all, HLbound, optRound_KO, Parallelisation, FinalFolderName);  
end

%%% Discriminate fast effects from long-term effects
X = estim.Results.KnockOutNodes.BIC_values;
Y = estim.Results.KnockOutNodesFast.BIC_values;
M1 = min(X,Y); M2 = min(X,Y);
figure, plot(X,Y,'.k', 'MarkerSize', 15), hold on
plot([M1, M2], [M1, M2], '-k');


% === End of the script === %
