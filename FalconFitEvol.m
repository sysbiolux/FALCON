function [estim] = FalconFitEvol(varargin)
% FalconFitEvol re-optimise the model and plot the evolution of fitting costs (3 rounds)
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

estim=varargin{1};
IC_Dist=varargin{2};
FinalFolderName=varargin{3};

h = waitbar(0,'Please wait...');

% Resimulation 3 times and collect fitting costs
global Costs
k=FalconIC(estim,IC_Dist);
waitbar(1/3,h,'Running Fitting Evolution Round 1/3 ...')
[~,~]=FalconObjFunPlus(estim,k);
Cost1=Costs;

Costs=[];
k=FalconIC(estim,IC_Dist);
waitbar(2/3,h,'Running Fitting Evolution Round 2/3 ...')
[~,~]=FalconObjFunPlus(estim,k);
Cost2=Costs;

Costs=[];
k=FalconIC(estim,IC_Dist);
waitbar(3/3,h,'Running Fitting Evolution Round 3/3 ...')
[~,~]=FalconObjFunPlus(estim,k);
Cost3=Costs;
close(h);

% Arrange the fitting costs as a matrix and plot
PlotCosts=NaN(3,max([length(Cost1),length(Cost2),length(Cost3)]));
PlotCosts(1,1:length(Cost1))=Cost1;
PlotCosts(2,1:length(Cost2))=Cost2;
PlotCosts(3,1:length(Cost3))=Cost3;

hconv=figure;
plot(1:max([length(Cost1),length(Cost2),length(Cost3)]),PlotCosts)
title('SSE convergence during optimization')
xlabel('iteration'); ylabel('SSE');
saveas(hconv,[FinalFolderName, filesep, 'OptimiserConvergence'],'tif')
saveas(hconv,[FinalFolderName, filesep, 'OptimiserConvergence'],'fig')

estim.Results.FitEvol.PlotCosts = PlotCosts;
estim.Results.FitEvol.Cost1 = Cost1;
estim.Results.FitEvol.Cost2 = Cost2;
estim.Results.FitEvol.Cost3 = Cost3;

end