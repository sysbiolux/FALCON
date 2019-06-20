function [estim] = FalconKONodes(varargin)
% FalconKO creates new models for each knocked-out node by creating a dummy
% inhibiting node. Then calculates the fitness of the model and compares
% all KO models with BIC.
% [estim] = FalconKONodes(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO,FinalFolderName)
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

%fetching values from arguments
estim = varargin{1};
fxt_all = varargin{2};
HLbound = varargin{3};
optRound_KO = varargin{4};
Parallelisation = varargin{5};
ToSave = 0;
if nargin > 5
    Folder = varargin{6};
    ToSave = 1;
end

estim_orig = estim;
Param_original = estim.param_vector;
Interactions_original = estim.Interactions;
bestcosts = fxt_all(:,1);
fval_collect = [];
xval_collect = [];

%% BIC calculation
N = numel(estim.Output)-sum(sum(isnan(estim.Output)));
Nodes = estim.state_names;
Nodes(estim.Input_idx(1, :)) = [];
pn = length(Nodes);
p = numel(Param_original);

BIC_complete = N*log(bestcosts) + (log(N))*p; %BIC for base model

PreviousOptions = estim.options;
SSthresh = estim.SSthresh;

%%% parameter perturbation and refitting
thisfig = figure; hold on
set(gca, 'TickLabelInterpreter', 'none')
suptitle('Virtual Node KO');
BICs = BIC_complete;
if size(BICs,1) < optRound_KO
    BICs = [BICs ; NaN((optRound_KO - size(BICs,1)), 1)];
end
wb = waitbar(0,'Please wait...');
for counter = 1:pn
    waitbar(counter/pn, wb, sprintf('Running Nodes Knock-out Round %d out of %d ...', counter, pn))
    
    Interactions = Interactions_original;
    estim = estim_orig;
    thisNode = Nodes(counter);
    
    Interactions = [Interactions; ['ix', 'X_INHIB_X', '-|', thisNode, '1', 'N', 'D']];
    estim.state_names = [estim.state_names, 'X_INHIB_X'];
    
    estim.Input_idx = [estim.Input_idx, ones(size(estim.Input_idx, 1), 1) .* length(estim.state_names)];
    estim.Input = [estim.Input, ones(size(estim.Input, 1), 1)];
    
    MeasFile = FalconData2File(estim);
    FalconInt2File(Interactions, 'KDN_TempFile.txt')
    
    estim = FalconMakeModel('KDN_TempFile.txt', MeasFile, HLbound);
    estim.options = PreviousOptions;
    estim.SSthresh = SSthresh;
    fval_all = [];
    xval_all = [];
    
    if Parallelisation
        parfor counterround = 1:optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            [~,fval] = FalconObjFun(estim, k); %objective function
            fval_all = [fval_all; fval];
            xval_all = [xval_all; xval];
        end
    else
        for counterround = 1:optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            [~,fval] = FalconObjFun(estim, k); %objective function
            fval_all = [fval_all; fval];
            xval_all = [xval_all; xval];
        end
    end
    fval_collect = [fval_collect, fval_all];
    xval_collect(:,:,counter) = xval_all;
    
    %% reduced model (- parameter)
    
    N_r = numel(estim.Output) - sum(sum(isnan(estim.Output))); %number of datapoints
    p_r = numel(estim.param_vector); %number of parameters
    %remove parameters related to the knocked-out node
    Is = Interactions_original(strcmp(Interactions_original(:, 2),thisNode), 5); %fetch outgoing parameters
    Is = [Is;Interactions_original(strcmp(Interactions_original(:, 4),thisNode), 5)]; %fetch incoming parameters
    for cc = length(Is):-1:1 %remove numbers
        if str2double(char(Is(cc))) >= 0
            Is(cc) = [];
        end
    end
    
    BICs = [BICs, N_r .* log(fval_all) + (log(N_r)) .* p_r];    
    %%Plot BIC values
    
    set(0, 'CurrentFigure', thisfig);
    figko = thisfig; hold on;
    
    boxplot(BICs), hold on
    sinaplot(BICs)
        
    set(gca, 'XTick', 1:length(Nodes)+1)
    Xtitles = ['Full'; Nodes'];
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
        saveas(figko, [Folder,filesep, 'KO_Nodes'], 'fig')        
    end
    
end
close(wb);

if ToSave
    saveas(figko, [Folder,filesep, 'KO_Nodes'], 'tif')
    saveas(figko, [Folder,filesep, 'KO_Nodes'], 'fig')
    saveas(figko, [Folder,filesep, 'KO_Nodes'], 'jpg')
    saveas(figko, [Folder,filesep, 'KO_Nodes'], 'svg')
    saveas(figko, [Folder,filesep, 'KO_Nodes'], 'pdf')
end

estim = estim_orig;

estim.Results.KnockOutNodes.Parameters = Xtitles';
estim.Results.KnockOutNodes.BIC_values = BICs;
estim.Results.KnockOutNodes.AllEvals = fval_collect;
estim.Results.KnockOut.AllValues = xval_collect;
delete('KDN_TempFile.txt')

end
