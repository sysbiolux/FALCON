function estim = FalconOptimize(estim)
% Wrap-up function to perform one round of optimizations, retrieving results and re-simulating the system %
% :: Input values ::
% estim             model variable produced by FalconMakeModel or FalconMakeGlobalModel functions
%
% :: Output values::
% Result            structure with optimization results

% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu, seb@delandtsheer.com

toc_all = []; x_all = []; fval_all = []; AICs = []; MSEs = []; BICs = []; NparamsL = [];

if estim.Parallelisation
    disp('parallelized run starting...')
    parfor counter = 1:estim.optRound %'parfor' will run the parallel computing toolbox
        tic, k = FalconIC(estim); %initial conditions
        [xval, fval, AIC, MSE, BIC, Nparams] = FalconObjFun(estim, k); %objective function
        toc_all = [toc_all; toc]; x_all = [x_all; xval]; fval_all = [fval_all; fval];
        AICs = [AICs; AIC]; MSEs = [MSEs; MSE]; BICs = [BICs; BIC]; NparamsL = [NparamsL; Nparams];   
    end

else
    disp('single-core run starting...')
    for counter = 1:estim.optRound %'parfor' will run the parallel computing toolbox
        tic, k = FalconIC(estim); %initial conditions
        [xval, fval, AIC, MSE, BIC, Nparams] = FalconObjFun(estim, k); %objective function
        toc_all = [toc_all; toc]; x_all = [x_all; xval]; fval_all = [fval_all; fval];
        AICs = [AICs; AIC]; MSEs = [MSEs; MSE]; BICs = [BICs; BIC]; NparamsL = [NparamsL; Nparams];  
    end
    close(h);
end

beep; pause(0.5); beep;


param_vector = estim.param_vector;
FinalFolderName = estim.FinalFolderName;

% Extract best parameter set and calculate means & SDs
bestx_index = fval_all == min(fval_all);
if sum(bestx_index)>1 %if two evaluations have the exact same cost we take one randomly
    bestx_index = bestx_index.*rand(size(bestx_index));
    bestx_index = bestx_index == max(bestx_index);
end
bestf = fval_all(bestx_index);
bestx = x_all(bestx_index, :);
bestt = toc_all(bestx_index);
bestAIC = AICs(bestx_index);
bestMSE = MSEs(bestx_index);
bestBIC = BICs(bestx_index);
bestNparams = NparamsL(bestx_index);
meanf = mean(fval_all);
stdf = std(fval_all, 0, 1);
meanx = mean(x_all, 1);
stdx = std(x_all, 0, 1);
meant = mean(toc_all);
stdt = std(toc_all, 0, 1);


% Display results in table format

% First table: Fitting Cost and Time

Heading = cell(1, 4);
Heading(1, 1) = {'-------'};
Heading(1, 2) = {'best'};
Heading(1, 3) = {'mean'};
Heading(1, 4) = {'S.D.'};
Heading_1stCol = {'MSE';'Time(s)'};
FitCost_Time = [Heading_1stCol num2cell([bestf bestt]') num2cell([meanf meant]') num2cell([stdf stdt]')];

disp('Summarized fitting cost and time:')
disp(' ')
disp([Heading; FitCost_Time])
disp(' ')

% Second Table: Parameter values

% Convert parameter names from symbolic to string
param_string = {};
for counter = 1:size(param_vector, 1)
    param_string(counter) = {char(param_vector(counter))};
end

Heading = cell(1,4);
Heading(1,1) = {'parameters'};
Heading(1,2) = {'best'};
Heading(1,3) = {'mean'};
Heading(1,4) = {'S.D.'};

param_name_best_mean_std = [param_string' num2cell(bestx') num2cell(meanx') num2cell(stdx')];

disp('Summarized identified parameters:')
disp(' ')
disp([Heading; param_name_best_mean_std])

if exist(FinalFolderName)
    useexcel = isExcelPresent();
    if useexcel        
        xlswrite([pwd filesep FinalFolderName filesep 'Summary_Optimised_Parameters.xls'], [Heading; param_name_best_mean_std])
    else
        tab = table(param_vector, bestx', meanx', stdx', 'VariableNames', {'Parameter', 'Best', 'Average', 'Std'});
        writetable(tab, [FinalFolderName, filesep, 'Summary_Optimised_Parameters.csv'], 'Delimiter', ',');
    end
end

%writing log file:
if exist(FinalFolderName)
    Start = FinalFolderName(9:end);
else
    Start = '?';
end

Now = datestr(now);
Now = strrep(Now, ':', '.');
fid = fopen([FinalFolderName, filesep, 'FALCONrun.log'], 'at');
fprintf(fid, 'FALCON log file \n');
fprintf(fid, 'Start timestamp: %s \n', Start);
fprintf(fid, 'Finish timestamp: %s \n', Now);
fprintf(fid, 'Final MSE: %d \n', bestf);
fprintf(fid, 'Network: \n');
[L,~] = size(estim.Interactions);
for c = 1:L
    fprintf(fid, '%s %s %s %s %s %s %s \n', estim.Interactions{c, :});
end

fprintf(fid, 'Dataset: \n');
[L,C] = size(estim.Input);
FormatIn = repmat('%d ', 1, C);
fprintf(fid, 'Inputs: \n');
for c = 1:L
    fprintf(fid, [FormatIn, '\n'], estim.Input(c, :));
end

[L,C] = size(estim.Input_idx);
FormatIn = repmat('%d ', 1, C);
fprintf(fid, 'Inputs indices: \n');
for c = 1:L
    fprintf(fid, [FormatIn, '\n'], estim.Input_idx(c, :));
end

[L,C] = size(estim.Output);
FormatIn = repmat('%d ', 1, C);
fprintf(fid, 'Outputs: \n');
for c = 1:L
    fprintf(fid, [FormatIn, '\n'], estim.Output(c,:));
end

[L,C] = size(estim.Output_idx);
FormatIn = repmat('%d ', 1, C);
fprintf(fid, 'Outputs indices: \n');
for c = 1:L
    fprintf(fid, [FormatIn, '\n'], estim.Output_idx(c,:));
end

fclose(fid);

disp('==============================================')
disp(' ')

estim.MaxTime = meant * 3;
estim.Results.Optimization.FittingCost = fval_all(:, 1);
estim.Results.Optimization.FittingTime = toc_all(:, end);
estim.Results.Optimization.ParamNames = estim.param_vector;
estim.Results.Optimization.BestParams = bestx;
estim.Results.Optimization.StateNames = estim.state_names;
estim.Results.Optimization.MSE = bestMSE;
estim.Results.Optimization.AIC = bestAIC;
estim.Results.Optimization.BIC = bestBIC;
estim.Results.Optimization.Nparams = bestNparams;

end















