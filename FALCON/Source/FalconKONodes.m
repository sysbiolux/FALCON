function [estim] = FalconKONodes(varargin)
% FalconKO creates new models for each knocked-out node by creating a dummy
% inhibiting node. Then calculates the fitness of the model and compares
% all KO models with BIC.
% [estim] = FalconKONodes(estim, fxt_all, HLbound, optRound_KO, Parallelisation, FinalFolderName)
%
% :: Input values ::
% estim             complete model definition
%
% :: Output value(s) ::
% estim             updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

%fetching values from arguments
estim = varargin{1};
optRound_KO = estim.optRound_KO;
Parallelisation = estim.Parallelisation;
Folder = estim.FinalFolderName;
ToSave = 1;

if ~isfield(estim.Results.Optimization, 'StateValueAll')
    estim = FalconSimul(estim, 0);
end
estim_orig = estim;
Param_original = estim.param_vector;
Interactions_original = estim.Interactions;
fval_collect = [];
nodeval_collect = [];
Diffs_collect = [];
%% BIC calculation
N = numel(estim.Output)-sum(sum(isnan(estim.Output)));
nOut = size(estim.Output, 2);
outNodes = estim.state_names(estim.Output_idx(1,:));
Nodes = estim.state_names;
Nodes(estim.Input_idx(1, :)) = [];
pn = length(Nodes);
p = numel(Param_original);

BIC_complete = N.*log(min(estim.Results.Optimization.FittingCost)) + (log(N)).*p; %BIC for base model
Diff_complete = mean(abs(estim.Results.Optimization.StateValueAll(:, estim.Output_idx(1,:)) - estim.Output));

%%% parameter perturbation and refitting
figko = figure; hold on;
set(gca, 'TickLabelInterpreter', 'none')
suptitle('Virtual Node KO');
BICs = BIC_complete;
SplitBICs = N.*log(Diff_complete) + (log(N)).*p;

for counter = 1:pn
    
    Interactions = Interactions_original;
    estim = estim_orig;
    thisNode = Nodes(counter);
    
    Interactions = [Interactions; ['ix', 'Xx_INHIB_xX', '-|', thisNode, '1', 'N', 'D']];
    estim.state_names = [estim.state_names, 'Xx_INHIB_xX'];
    
    estim.Input_idx = [estim.Input_idx, ones(size(estim.Input_idx, 1), 1) .* length(estim.state_names)];
    estim.Input = [estim.Input, ones(size(estim.Input, 1), 1)];
    
    for out = 1:nOut %if tested node is measured, we remove it from the measurements
        if strcmp(outNodes(out), thisNode)
            estim.Output(:,out) = nan(size(estim.Output, 2), 1);
        end
    end
    
    MeasFile = FalconData2File(estim);
    FalconInt2File(Interactions, 'KDN_TempFile.txt')
    
    estim = FalconMakeModel('KDN_TempFile.txt', MeasFile);
    estim.options = estim_orig.options;
    estim.SSthresh = estim_orig.SSthresh; estim.ObjFunction = estim_orig.ObjFunction;
    estim.Results.Optimization.BestParams = estim_orig.Results.Optimization.BestParams;
    fxt_all = [];
    
    if Parallelisation
        parfor counterround = 1:optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            tic
            [xval,fval] = FalconObjFun(estim, k); %objective function
            fxt_all = [fxt_all; [fval xval toc]];
        end
    else
        for counterround = 1:optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            tic
            [xval,fval] = FalconObjFun(estim, k); %objective function
            fxt_all = [fxt_all; [fval xval toc]];
        end
    end
    minf = find(fxt_all(:,1)==min(fxt_all(:,1)));
    minf = minf(1);
    bestx = fxt_all(minf, 2:end-1);
    bestf = fxt_all(minf, 1);

    [~, nodevals, fval] = FalconSimul(estim, 0);
    fval_collect = [fval_collect, fval];
    nodeval_collect = [nodeval_collect; nodevals];
    InNodes = nodevals(:, estim.Output_idx(1,:));
    Diffs = mean(abs(InNodes - estim.Output));
    Diffs_collect = [Diffs_collect; Diffs];
    
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
    
    BICs = [BICs, N_r .* log(bestf) + (log(N_r)) .* p_r];
    SplitBICs = [SplitBICs; N_r .* log(Diffs) + (log(N_r)) .* p_r];
    
    %%Plot BIC values
    
    figure(figko)
    isGreen = min(BICs) < BICs(1);
    subplot(2,1,1)
    b = bar(min(BICs, [], 1), 'r'); hold on
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.5 0.5 0.5];
    IdxGreen = find(isGreen);
    for Green = 1:length(IdxGreen)
        b.CData(IdxGreen(Green),:) = [0 1 0];
    end
    hold on,
    set(gca, 'XTick', 1:length(Nodes)+1)
    Xtitles = ['Full'; Nodes'];
    set(gca, 'XTicklabel', Xtitles);
    xlabel('');
    set(gca, 'XTickLabelRotation', 45)
    ylabel('BIC');
    hold off
    Min = min(BICs(:)); 
    Max = max(BICs(:));
    axis([0.5 counter+1.5 Min-0.1*abs(Min) Max+0.1*abs(Max)])
    drawnow;
    subplot(2,1,2)
    imagesc(SplitBICs');
    set(gca, 'ytick', 1:size(SplitBICs,2))
    set(gca, 'xtick', 1:(pn+1))
    set(gca, 'yticklabel',estim.state_names(estim.Output_idx(1, :)))
    set(gca, 'xticklabel',['Full', Nodes])
    set(gca, 'XTickLabelRotation', 45)

    colorbar, drawnow;
    
    if ToSave
        saveas(figko, [Folder,filesep, 'KO_Nodes'], 'fig')        
    end
    
end

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
estim.Results.KnockOutNodes.AllValues = nodeval_collect;
estim.Results.KnockOutNodes.AllDiffs = Diffs_collect;
estim.Results.KnockOutNodes.SplitBICs = SplitBICs;
delete('KDN_TempFile.txt')

end
