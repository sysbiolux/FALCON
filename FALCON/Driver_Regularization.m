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
InputFile = 'ToyDiff.xlsx'; %for xls
GlobalEdgesList = 'ToyDiff_fixed.xlsx';
ContextsList = {'CL1a','CL2a','CL3a','CL4a'};
MeasFileList = {};
for f = 1:length(ContextsList)
    MeasFileList = [MeasFileList, ['ToyDiff_meas_', char(ContextsList(f)), '.xlsx']];
end

% Builds a FALCON model for optimisation
estim = FalconMakeGlobalModel(InputFile, GlobalEdgesList, MeasFileList, ContextsList); %make the model

%%%%
% modify optimization parameters here
estim = FalconPreProcess(estim, 'noise', 5);

%%%%

% Optimization
estim = FalconOptimize(estim);

%%% Re-simulate results based on the best optimised parameter set
[estim, StateValues, MSE] = FalconSimul(estim);


%%% Select Regularization
% estim.Reg = 'L1'; % sum of parameters
% estim.Reg = 'L1Groups'; % sum of intra-group distances
estim.Reg = 'LCluster'; % spatial non-uniformity
% estim.Reg = 'L1Smooth'; % sequential distances
% estim.Reg = 'L2'; % sum of squared parameters
% estim.Reg = 'L1/2'; % sum of parameters square roots
% estim.Reg = 'Ldrug'; % L1/2 and L1Groups
% estim.Reg = 'PruneCluster'; % L1/2 and LCluster
% estim.Reg = 'LTriple'; % L1/2 and LCluster and L1Smooth


%%

% fix Lambda values here!
PowerLambda1 = -12:0.5:-2; % (L1/2 in combined schemes)
PowerLambda2 = -12:1:-2; % (L1group or LCluster in combined schemes)
PowerLambda3 = -2:2:2; % (L1Smooth when present)


ListLambda1 = [0,2.^PowerLambda1];
ListLambda2 = [0,2.^PowerLambda2];
ListLambda3 = [0,2.^PowerLambda3];

estim.RegMatrix.Groups = (reshape(1:length(estim.param_vector), length(ContextsList), length(estim.param_vector)/length(ContextsList)))';
estim.RegMatrix.Cluster = estim.RegMatrix.Groups;
estim.RegMatrix.Smooth = estim.RegMatrix.Groups;

Results = [];

%% this for 1-dimensional regularization


for L = ListLambda1
    estim.Lambda = L;
    % Optimization
    estim = FalconOptimize(estim);
    Results = [Results; ...
        min(estim.Results.Optimization.FittingCost), ...
        estim.Results.Optimization.BestParams, ...
        min(estim.Results.Optimization.FittingTime), ...
        estim.Results.Optimization.MSE, ...
        estim.Results.Optimization.BIC, ...
        estim.Results.Optimization.AIC, ...
        estim.Results.Optimization.Nparams];
end

L = length(ListLambda1);
figure,
subplot(3,2,1)
plot(Results(:,2:end-5)- repmat(Results(1,2:end-5),size(Results,1),1))
set(gca,'XTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
ylabel('bias')
subplot(3,2,2)
plot(Results(:,2:end-5))
set(gca,'XTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
ylabel('k')
subplot(3,2,3)
plot(1:L, log10(Results(:,end-3)), 'b'), hold on
plot(2:L, log10(Results(2:end,1) - Results(2:end,end-3)), 'm'), hold on
plot(1:L, log10(Results(:,1)), 'k'), hold on
set(gca,'XTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
ylabel('costs')
legend({'log(MSE)', 'log(Reg)', 'log(Total)'})
subplot(3,2,4)
plot(1:L, Results(:,end-2)), hold on
idx = find(Results(:,end-2) == min(Results(:,end-2)));
plot(idx, Results(idx,end-2), '*r')
set(gca,'XTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
legend({'BIC'})
subplot(3,2,5)
plot(1:L, Results(:,end-1))
set(gca,'XTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
legend({'AIC'})
subplot(3,2,6)
plot(1:L, Results(:,end))
set(gca,'XTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
legend({'#params'})



%% this for 2-dimensional regularization
estim.Reg = 'PruneCluster';
Results = [];
PowerLambda1 = -12:1:-2; % (L1/2 in combined schemes)
PowerLambda2 = -8:1:-2; % (L1group or LCluster in combined schemes)

ListLambda1 = [0,2.^PowerLambda1];
ListLambda2 = [0,2.^PowerLambda2];

for L1 = ListLambda1
    for L2 = ListLambda2
        estim.Lambda = [L1, L2];
        % Optimization
        estim = FalconOptimize(estim);
        Results = [Results; ...
            min(estim.Results.Optimization.FittingCost), ...
            estim.Results.Optimization.BestParams, ...
            min(estim.Results.Optimization.FittingTime), ...
            estim.Results.Optimization.MSE, ...
            estim.Results.Optimization.BIC, ...
            estim.Results.Optimization.AIC, ...
            estim.Results.Optimization.Nparams];
    end
end

L = length(ListLambda1);

BICs = reshape(Results(:,end-2), length(ListLambda1), length(ListLambda2));
MSEs = reshape(Results(:,end-3), length(ListLambda1), length(ListLambda2));
AICs = reshape(Results(:,end-1), length(ListLambda1), length(ListLambda2));
NParams = reshape(Results(:,end), length(ListLambda1), length(ListLambda2));

figure,
subplot(2,2,1)
imagesc(BICs)
set(gca,'XTick',1:length(ListLambda2))
set(gca,'YTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
set(gca,'yticklabels',{'base'; num2str(PowerLambda1(:))})
title('BIC')
subplot(2,2,2)
imagesc(AICs)
set(gca,'XTick',1:length(ListLambda2))
set(gca,'YTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
set(gca,'yticklabels',{'base'; num2str(PowerLambda1(:))})
title('AIC')
subplot(2,2,3)
imagesc(MSEs)
set(gca,'XTick',1:length(ListLambda2))
set(gca,'YTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
set(gca,'yticklabels',{'base'; num2str(PowerLambda1(:))})
title('MSE')
subplot(2,2,4)
imagesc(NParams)
set(gca,'XTick',1:length(ListLambda2))
set(gca,'YTick',1:length(ListLambda1))
set(gca,'xticklabels',{'base'; num2str(PowerLambda1(:))})
set(gca,'xticklabelrotation',90)
set(gca,'yticklabels',{'base'; num2str(PowerLambda1(:))})
title('#params')



