function [estim] = FalconFitEvol(varargin)
% FalconFitEvol re-optimise the model and plot the evolution of fitting costs (3 rounds)
% Cannot run in parallel because we need a global variable, as we ask
% fmincon to return additional info.
% FalconFitEvol(estim,IC_Dist,[FinalFolderName])
% 
% :: Input values ::
% estim             complete model definition
% IC_Dist           type of distribution for the initial condition 
% FinalFolderName   name of the folder for saving results (optional)
%
% :: Output values ::
% estim             updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu
estim = varargin{1};
FinalFolderName = estim.FinalFolderName;

% Resimulation 3 times and collect fitting costs
global Costs
Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost1 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost2 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost3 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost4 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost5 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost6 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost7 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost8 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost9 = Costs;

Costs = [];
k = FalconIC(estim);
[~,~] = FalconObjFunPlus(estim, k);
Cost10 = Costs;


% Arrange the fitting costs as a matrix and plot
PlotCosts = NaN(10, max([length(Cost1), length(Cost2), length(Cost3),length(Cost4), length(Cost5), length(Cost6)...
    length(Cost7), length(Cost8), length(Cost9), length(Cost10)]));
PlotCosts(1, 1:length(Cost1)) = Cost1;
PlotCosts(2, 1:length(Cost2)) = Cost2;
PlotCosts(3, 1:length(Cost3)) = Cost3;
PlotCosts(4, 1:length(Cost4)) = Cost4;
PlotCosts(5, 1:length(Cost5)) = Cost5;
PlotCosts(6, 1:length(Cost6)) = Cost6;
PlotCosts(7, 1:length(Cost7)) = Cost7;
PlotCosts(8, 1:length(Cost8)) = Cost8;
PlotCosts(9, 1:length(Cost9)) = Cost9;
PlotCosts(10, 1:length(Cost10)) = Cost10;

PlotCosts = log10(PlotCosts);

hconv = figure;
plot(1:max([length(Cost1), length(Cost2), length(Cost3),length(Cost4), length(Cost5), length(Cost6)...
    length(Cost7), length(Cost8), length(Cost9), length(Cost10)]), PlotCosts)
title('Convergence during optimization')
xlabel('iteration'); ylabel('log_{10}(MSE)');
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'tif')
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'fig')
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'jpg')
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'svg')

estim.Results.FitEvol.PlotCosts = PlotCosts;

end