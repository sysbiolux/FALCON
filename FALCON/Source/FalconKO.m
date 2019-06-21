function [estim] = FalconKO(varargin)
% FalconKO creates new models for each knock-outed parameter by setting parameter value to 0 and re-optimise
% [estim] = FalconKO(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO,FinalFolderName)
%
% :: Input values ::
% estim             complete model definition
% fxt_all           all fitting costs, parameter values and running time during the different optimisations
% MeasFile          experimental data
% HLbound           qualitative threshold between high and low range of parameter values
% optRound_KO       number of optimization rounds
% Parallelisation   whether parallel computing is used
% FinalFolderName   name of the folder for saving results
%
% :: Output value(s) ::
% estim             updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fetching values from arguments
estim = varargin{1};
fxt_all = varargin{2};
HLbound = varargin{3};
optRound_KO = varargin{4};
Parallelisation = varargin{5};
ToSave = 0;
if nargin > 5
    FinalFolderName = varargin{6};
    ToSave = 1;
end

estim_orig = estim;
MeasFile = FalconData2File(estim);
Param_original = estim.param_vector;
Interactions_original = estim.Interactions;
num_params = estim.param_vector;
bestcosts = fxt_all(:, 1);
fval_collect = [];
nodeval_collect = [];

%% BIC calculation
N = numel(estim.Output) - sum(sum(isnan(estim.Output)));
p= numel(Param_original);

BIC_complete = N*log(bestcosts) + (log(N))*p; %BIC for base model

PreviousOptions = estim.options;
SSthresh = estim.SSthresh;

%%% parameter perturbation and refitting
thisfig = figure; hold on
set(gca, 'TickLabelInterpreter', 'none')
suptitle('Virtual Interaction KO');
BICs = BIC_complete;
if size(BICs,1) < optRound_KO
    BICs = [BICs ; NaN((optRound_KO - size(BICs,1)), 1)];
end
wb = waitbar(0, 'Please wait...');
for counter = 1:p
    waitbar(counter/p, wb, sprintf('Running Knock-out Round %d out of %d ...', counter, p))
    
    Interactions = Interactions_original;
    replace_idx = find(ismember(Interactions_original(:, 5),Param_original(counter))); %index of interactions to modify
    disp('Removing interaction...')
    disp(Interactions_original(replace_idx, :))
    for counter2 = 1:length(replace_idx) %for each one
        Interactions{replace_idx(counter2), 5} = '0'; %replace weight by '0'
    end
    FalconInt2File(Interactions, 'KD_TempFile.txt') %write this as a temp file
    estim = FalconMakeModel('KD_TempFile.txt', MeasFile, HLbound); %make model variant
    estim.options = PreviousOptions; %optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',3000,'MaxIter',3000); % Default
    estim.SSthresh = SSthresh;
    fval_all = [];
    nodeval_all = [];
    
    if Parallelisation
        parfor counterround = 1:optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            [xval,fval] = FalconObjFun(estim, k); %objective function
            fval_all = [fval_all; fval];
            [nodevals, ~, ~, ~, ~] = FalconSimul(estim, xval, [0 0 0 0 0]);
            nodeval_all = [nodeval_all; nodevals];
        end
    else
        for counterround = 1:optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            [xval,fval] = FalconObjFun(estim, k); %objective function
            fval_all = [fval_all; fval];
            [nodevals, ~, ~, ~, ~] = FalconSimul(estim, xval, [0 0 0 0 0]);
            nodeval_all = [nodeval_all; nodevals];
        end
    end
    fval_collect = [fval_collect, fval_all];
    nodeval_collect = [nodeval_collect; nodeval_all];
    
    %% reduced model (-1 parameter)
    
    N_r = numel(estim.Output) - sum(sum(isnan(estim.Output))); %number of datapoints
    p_r = (numel(estim.param_vector)); %number of parameters
    
    BICs = [BICs, N_r .* log(fval_all) + (log(N_r)) .* p_r];
    
    %%Plot BIC values
    
    figko = thisfig;
    set(0, 'CurrentFigure', thisfig); %hold on;
    
    boxplot(BICs), hold on
    sinaplot(BICs)
        
    set(gca, 'XTick', 1:(length(num_params)+1))
    Xtitles = ['Full'; Param_original];
    set(gca, 'XTicklabel', Xtitles);

    xlabel('');
    set(gca, 'XTickLabelRotation', 45)
    ylabel('Bayesian Information Criterion (BIC)');
    hold off
    Min = min(BICs(:));
    Max = max(BICs(:));
    axis([0.5 counter+1.5 Min-0.1*abs(Min) Max+0.1*abs(Max)])
    drawnow;
    if ToSave
        saveas(figko, [FinalFolderName,filesep, 'KO_interactions'], 'fig')
    end
    
end
close(wb);

if ToSave
    saveas(figko, [FinalFolderName,filesep, 'KO_interactions'], 'tif')
    saveas(figko, [FinalFolderName,filesep, 'KO_interactions'], 'fig')
    saveas(figko, [FinalFolderName,filesep, 'KO_interactions'], 'jpg')
    saveas(figko, [FinalFolderName,filesep, 'KO_interactions'], 'svg')
    saveas(figko, [FinalFolderName,filesep, 'KO_interactions'], 'pdf')
end

estim = estim_orig;

estim.Results.KnockOut.Parameters = Xtitles';
estim.Results.KnockOut.BIC_values = BICs;
estim.Results.KnockOut.AllEvals = fval_collect;
estim.Results.KnockOut.AllValues = nodeval_collect;
delete('KD_TempFile.txt')

end
