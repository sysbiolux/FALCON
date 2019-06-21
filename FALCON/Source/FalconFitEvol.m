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
IC_Dist = varargin{2};
FinalFolderName = varargin{3};

h = waitbar(0,'Please wait...');

% Resimulation 3 times and collect fitting costs
global Costs
Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(1/10, h, 'Running Fitting Evolution Round 1/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost1 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(2/10, h, 'Running Fitting Evolution Round 2/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost2 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(3/10, h, 'Running Fitting Evolution Round 3/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost3 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(4/10, h, 'Running Fitting Evolution Round 4/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost4 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(5/10, h, 'Running Fitting Evolution Round 5/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost5 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(6/10, h, 'Running Fitting Evolution Round 6/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost6 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(7/10, h, 'Running Fitting Evolution Round 7/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost7 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(8/10, h, 'Running Fitting Evolution Round 8/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost8 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(9/10, h, 'Running Fitting Evolution Round 9/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost9 = Costs;

Costs = [];
k = FalconIC(estim, IC_Dist);
waitbar(10/10, h, 'Running Fitting Evolution Round 10/10 ...')
[~,~] = FalconObjFunPlus(estim, k);
Cost10 = Costs;

close(h);

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

hconv = figure;
plot(1:max([length(Cost1), length(Cost2), length(Cost3),length(Cost4), length(Cost5), length(Cost6)...
    length(Cost7), length(Cost8), length(Cost9), length(Cost10)]), PlotCosts)
title('Convergence during optimization')
xlabel('iteration'); ylabel('MSE');
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'tif')
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'fig')
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'jpg')
saveas(hconv, [FinalFolderName, filesep, 'OptimiserConvergence'],'svg')

estim.Results.FitEvol.PlotCosts = PlotCosts;

end