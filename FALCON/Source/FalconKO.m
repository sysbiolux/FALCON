function [estim] = FalconKO(varargin)
% FalconKO creates new models for each knock-outed parameter by setting parameter value to 0 and re-optimise
% [estim] = FalconKO(estim)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fetching values from arguments
estim = varargin{1};

estim_orig = estim;
MeasFile = FalconData2File(estim);
Param_original = estim.param_vector;
Interactions_original = estim.Interactions;
num_params = estim.param_vector;
fval_collect = [];
nodeval_collect = [];

%% BIC calculation
N = numel(estim.Output) - sum(sum(isnan(estim.Output)));
p= numel(Param_original);

BIC_complete = N.*log(min(estim.Results.Optimization.FittingCost)) + (log(N)).*p; %BIC for base model

PreviousOptions = estim.options;
SSthresh = estim.SSthresh;

%%% parameter perturbation and refitting
thisfig = figure; hold on
set(gca, 'TickLabelInterpreter', 'none')
sgtitle('Virtual Interaction KO');
BICs = BIC_complete;
if size(BICs,1) < estim.optRound_KO
    BICs = [BICs ; NaN((estim.optRound_KO - size(BICs,1)), 1)];
end

for counter = 1:p
    Interactions = Interactions_original;
    replace_idx = find(ismember(Interactions_original(:, 5),Param_original(counter))); %index of interactions to modify
    disp('Removing interaction...')
    disp(Interactions_original(replace_idx, :))
    for counter2 = 1:length(replace_idx) %for each one
        Interactions{replace_idx(counter2), 5} = '0'; %replace weight by '0'
    end
    FalconInt2File(Interactions, 'KD_TempFile.txt') %write this as a temp file
    estim = FalconMakeModel('KD_TempFile.txt', MeasFile); %make model variant
    estim.options = PreviousOptions; %optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',3000,'MaxIter',3000); % Default
    estim.SSthresh = SSthresh; estim.ObjFunction = estim_orig.ObjFunction;
    b = estim_orig.Results.Optimization.BestParams;
    if counter == 1
        b = b(2:end);
    elseif counter == p
        b = b(1:end-1);
    else
        b = [b(1:(counter-1)), b((counter+1):end)];
    end
    estim.Results.Optimization.BestParams = b;
    if isfield(estim, 'Weights')
        estim.Weights = estim_orig.Weights;
    end
    fval_all = [];
    nodeval_all = [];
    
    if estim.Parallelisation
        parfor counterround = 1:estim.optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            [xval,fval] = FalconObjFun(estim, k); %objective function
            fval_all = [fval_all; fval];
            [~, nodevals, ~] = FalconSimul(estim, 0);
            nodeval_all = [nodeval_all; nodevals];
        end
    else
        for counterround = 1:estim.optRound_KO %computation for modified model
            k = FalconIC(estim); %initial conditions
            [xval,fval] = FalconObjFun(estim, k); %objective function
            fval_all = [fval_all; fval];
            [~, nodevals, ~] = FalconSimul(estim, 0);
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
    isGreen = min(BICs)<min(BICs(:,1));
    
    b = bar(min(BICs), 'r'); hold on
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.5 0.5 0.5];
    IdxGreen = find(isGreen);
    for Green = 1:length(IdxGreen)
        b.CData(IdxGreen(Green),:) = [0 1 0];
    end
    hold on,
    sinaplot(BICs)
        
    set(gca, 'XTick', 1:(length(num_params)+1))
    Xtitles = ['Full'; Param_original];
    set(gca, 'XTicklabel', Xtitles);

    xlabel('');
    set(gca, 'XTickLabelRotation', 45)
    ylabel('BIC');
    hold off
    Min = min(BICs(:));
    Max = max(BICs(:));
    axis([0.5 counter+1.5 Min-0.1*abs(Min) Max+0.1*abs(Max)])
    drawnow;

    saveas(figko, [estim_orig.FinalFolderName,filesep, 'KO_interactions'], 'fig') %save temporary figure in case of crash


end

estim = estim_orig;

saveas(figko, [estim.FinalFolderName,filesep, 'KO_interactions'], 'tif')
saveas(figko, [estim.FinalFolderName,filesep, 'KO_interactions'], 'fig')
saveas(figko, [estim.FinalFolderName,filesep, 'KO_interactions'], 'jpg')
saveas(figko, [estim.FinalFolderName,filesep, 'KO_interactions'], 'svg')
saveas(figko, [estim.FinalFolderName,filesep, 'KO_interactions'], 'pdf')



estim.Results.KnockOut.Parameters = Xtitles';
estim.Results.KnockOut.BIC_values = BICs;
estim.Results.KnockOut.AllEvals = fval_collect;
estim.Results.KnockOut.AllValues = nodeval_collect;
delete('KD_TempFile.txt')

end
