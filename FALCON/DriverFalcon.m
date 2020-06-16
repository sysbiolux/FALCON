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
%                   Current version: v1.3                      %
%                      (July 2020)                             %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all % clear screen and workspace 

% Model and measurement files
InputFile = 'PDGF_model.xlsx'; %for xls
MeasFile = 'PDGF_meas.xlsx'; %for xls

% InputFile='PDGF_model.sif'; %for sif
% MeasFile={'PDGF_meas_in.csv';'PDGF_meas_out.csv';'PDGF_meas_err.csv'}; %for csv

% Builds a FALCON model for optimisation
estim = FalconMakeModel(InputFile, MeasFile); %make the model

%%%%
% modify optimization parameters here
%%%%

% Optimization
estim = FalconOptimize(estim);

%%% Re-simulate results based on the best optimised parameter set
[estim, StateValues, MSE] = FalconSimul(estim);

save([estim.FinalFolderName filesep 'OptimizedModel'])

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
    estim = FalconResample(estim);
end

%%% Sensitivity analysis
if estim.LPSA_Analysis 
    estim = FalconLPSA(estim);
end

%%% Knock-out analysis
if estim.KO_Analysis
    estim = FalconKO(estim);
end

%%% Nodes Knock-out analysis
if estim.KO_Nodes_Analysis
    estim = FalconKONodes(estim);
end

%%% Fast KO
if estim.KO_Nodes_fast
    estim = FalconKONodes_fast(estim);  
end
% 
% %%% Nodes Knock-out efficency analysis
% if estim.KO_Nodes_Analysis_eff
%     estim=FalconKONodes_eff(estim, fxt_all, efficiency_range, HLbound, optRound_KO, FinalFolderName);    
% 
% end

%%% Discriminate fast effects from long-term effects
X = estim.Results.KnockOutNodes.BIC_values;
Y = estim.Results.KnockOutNodesFast.BIC_values;
figure, plot(X,Y,'.k', 'MarkerSize', 15), hold on
M1 = min(X,Y); M2 = min(X,Y);
plot([M1, M2], [M1, M2], '-k');
Nodes = estim.state_names;
Nodes(estim.Input_idx(1, :)) = [];
text(X(1), Y(1), 'base')
for tx = 1:length(X)-1
    text(X(tx+1), Y(tx+1), Nodes(tx))
end
xlabel('BIC for re-optimized model')
ylabel('BIC for non-adapted model')




% === End of the script === %
