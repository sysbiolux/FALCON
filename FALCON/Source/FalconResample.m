function [estim, Ks, Costs] = FalconResample(varargin)
% FalconResample creates a new dataset based on the mean and standard deviation of the original dataset.
% Useful to assess the variability of the parameter estimates, in conjunction with FalconLPSA
% [Ks,Costs,estim]=FalconResample(estim, bestx, optRound_Sampling, New_Data_Sampling, CV_cutoff, FinalFolderName)
%
% :: Input values ::
% estim               complete model definition
% bestx               the vector of optimal parameters
% optRound_Sampling   the number of times the resampling of experimental replicates is performed
% New_Data_Sampling   the number of experimental replicates from the original data
% CV_cutoff           the coefficient of variation (CV) cutoff value to determine if the optimised parameters have a large distribution
% FinalFolderName     name of the folder for saving results
%
% :: Output values ::
% Ks                  optimised parameter values from the optmisation based on the resampled data
% Costs               optimised fitting costs from the optmisation based on the resampled data
% estim               updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

estim = varargin{1};
bestx = estim.Results.Optimization.BestParams; 
N = estim.NDatasets;

estim_orig = estim;
Folder = estim.FinalFolderName;

Ks = [];
Costs = [];
Data = estim.Output;
SD = estim.SD;

parfor n = 1:N
    estim2 = estim;
    new = normrnd(Data,SD);
    % re-normalize
    maxnew = max(new(:));
    minnew = min(new(:));
    new = (new - minnew) ./ (maxnew - minnew);
    % new synthetic outputs
    estim2.Output = new;
    k = FalconIC(estim);
    [xval,fval] = FalconObjFun(estim2, k);
    Ks = [Ks; xval];
    Costs = [Costs; fval];
end

estim = estim_orig;


disp(['Cost for perturbated measurements: ', num2str(mean(Costs)), ' +- ', num2str(std(Costs))])
h1 = figure; hold on;
plot(bestx, '*r', 'MarkerSize', 12), hold on;
errorbar(mean(Ks), std(Ks), '.k', 'LineWidth', 1); hold on;
ylim([0 1])
title('Distribution of optimised parameters after resampling')
hold on;
set(gca, 'Xtick', 1:length(mean(Ks)), 'XTickLabel', estim.param_vector, 'XGrid', 'on');
set(gca, 'Xticklabelrotation', 90)
legend('Mean','Std')


saveas(h1, [Folder, filesep, 'Resampling'], 'tif');
saveas(h1, [Folder, filesep, 'Resampling'], 'fig');
saveas(h1, [Folder, filesep, 'Resampling'], 'jpg');
saveas(h1, [Folder, filesep, 'Resampling'], 'svg');


h1 = figure;
title('Distribution of optimised parameters after resampling')
hold on;
set(gca, 'Xtick', 1:length(mean(Ks)), 'XTickLabel', estim.param_vector, 'XGrid', 'on'); hold on
set(gca, 'Xticklabelrotation', 90), hold on
sinaplot(Ks);

saveas(h1, [Folder, filesep, 'ResamplingSinaplot'], 'tif');
saveas(h1, [Folder, filesep, 'ResamplingSinaplot'], 'fig');
saveas(h1, [Folder, filesep, 'ResamplingSinaplot'], 'jpg');
saveas(h1, [Folder, filesep, 'ResamplingSinaplot'], 'svg');


estim.Results.Resampling.Parameters = estim.param_vector';
estim.Results.Resampling.OptimisedParameter = Ks;
estim.Results.Resampling.OptimisedSD = std(Ks);
estim.Results.Resampling.Costs = Costs;

Heading = cell(1,3);
Heading(1,1) = {'parameters'};
Heading(1,2) = {'mean'};
Heading(1,3) = {'S.D.'};

Resampling = [estim.param_vector num2cell(mean(Ks, 1)') num2cell(std(Ks, 0, 1)')];

disp('Summary of resampling:')
disp(' ')
disp([Heading; Resampling])
disp(' ')


excelpresent = isExcelPresent();
if excelpresent
    xlswrite([Folder filesep 'Summary_Resampling.xls'], [Heading; Resampling]);
else
    setupxlwrite()
    xlwrite([Folder filesep 'Summary_Resampling.xls'], [Heading;Resampling]);
end


end