function [estim, cost_SA]=FalconLPSA(varargin)
% FalconLPSA performs a local parameters sensitivity analysis (LPSA) on an optimized model
% [estim, cost_SA]=FalconLPSA(estim, bestx, MeasFile, HLbound, optRound_LPSA, increment, option, Parallisation, FinalFolderName);
%
% :: Input values ::
% estim             complete model definition
%
% :: Output values ::
% estim             updated model definition
% cost_SA           fitting cost for each perturbation
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

estim = varargin{1};

estim_orig = estim;

%%%1) Resampling to choose cut-off
if estim.Fast_Option
    [~, ~, Costs] = FalconResample(estim);
    CutOff = mean(Costs) + std(Costs);
else
    CutOff = min(estim.Results.Optimization.FittingCost)*2;
end

p = estim.Results.Optimization.BestParams;
p_SA = zeros(2*estim.LPSA_Increments+1, length(p));
p_SA(estim.LPSA_Increments+1,:) = p;

for counter = 1:length(p)
    Min = estim.LB(counter) + (5 * eps);
    Max = estim.UB(counter) - (5 * eps);
    
    if ~isempty(estim.A)
        if sum(estim.A(:, counter)) > 0
            IdxC = find(estim.A(:, counter));
            Max = min(estim.b(IdxC), Max);
        end
    end
    
    for counter2 = 1:(estim.LPSA_Increments-1)
        p_pert_values_pos = p(counter) + ((Max-p(counter)) / estim.LPSA_Increments * counter2);
        p_pert_values_neg = p(counter) -((p(counter)-Min) / estim.LPSA_Increments * counter2);
        p_SA(estim.LPSA_Increments+1+counter2,counter) = p_pert_values_pos;
        p_SA(estim.LPSA_Increments+1-counter2,counter) = p_pert_values_neg;
    end
    
    p_SA(end,counter) = Max;
end

%%% Pick index and evaluate
cost_SA = NaN(size(p_SA));
cost_Error_SA = NaN(size(p_SA));
params_SA = NaN([size(p_SA), length(p)]);

%% Resimulation with perturbed parameter values
estim = estim_orig;
thisfig = figure; hold on
set(gca, 'TickLabelInterpreter', 'none')
Param_original = estim.param_vector;
Interactions_original = estim.Interactions;
MidPoint = ceil(size(p_SA, 1) / 2);

num_plots = length(Param_original);
NLines = ceil(sqrt(num_plots));
NCols = ceil(num_plots/NLines);

Ident_All = [];


for counter =  1:size(p_SA, 2) %for each parameter
    
    OK = 1;
    %1) for middle to 0
    for counter2 = MidPoint:-1:1
        if OK == 1
            Interactions = Interactions_original;
            replace_idx = find(ismember(Interactions_original(:, 5), Param_original(counter)));
            for counter3 = 1:length(replace_idx)
                Interactions{replace_idx(counter3), 5} = char(num2str(p_SA(counter2, counter)));
            end
            FalconInt2File(Interactions, 'LPSA_TempFile.txt');
            estim = FalconMakeModel('LPSA_TempFile.txt', estim.MeasFile);
            
            x_all = [];
            fval_all = [];
            
            if estim.Parallelisation
                parfor Round = 1:estim.optRound_LPSA %'parfor' will run the parallel computing toolbox
                    k = FalconIC(estim); %initial conditions
                    [xval,fval] = FalconObjFun(estim, k); %objective function
                    %collecting the run's outcome
                    x_all = [x_all; xval];
                    fval_all = [fval_all; fval];
                end
                
            else
                for Round = 1:estim.optRound_LPSA %'parfor' will run the parallel computing toolbox
                    k = FalconIC(estim); %initial conditions
                    [xval,fval] = FalconObjFun(estim, k); %objective function
                    %collecting the run's outcome
                    x_all = [x_all; xval];
                    fval_all = [fval_all; fval];
                end
            end
            
            cost_SA(counter2,counter) = min(fval_all);
            cost_Error_SA(counter2,counter) = std(fval_all);
            params_SA(counter2, counter, find(ismember(Param_original, estim.param_vector))) = ...
                x_all(min(find(ismember(fval_all, min(fval_all)))), :);
            if estim.Fast_Option
                if min(fval_all) > CutOff
                    OK = 0;
                end
            end
        end
        
    end
    
    %2) from middle to 1
    OK = 1;
    for counter2 =  MidPoint+1:size(p_SA, 1)
        
        if OK == 1
            Interactions = Interactions_original;
            replace_idx = find(ismember(Interactions_original(:, 5), Param_original(counter)));
            for counter3 = 1:length(replace_idx)
                Interactions{replace_idx(counter3), 5} = char(num2str(p_SA(counter2, counter)));
            end
            FalconInt2File(Interactions, 'LPSA_TempFile.txt');
            estim=FalconMakeModel('LPSA_TempFile.txt', estim.MeasFile);
            
            x_all = [];
            fval_all = [];
            
            if estim.Parallelisation == 1
                parfor Round = 1:estim.optRound_LPSA %'parfor' will run the parallel computing toolbox
                    k = FalconIC(estim); %initial conditions
                    [xval,fval] = FalconObjFun(estim, k); %objective function
                    %collecting the run's outcome
                    x_all = [x_all; xval];
                    fval_all = [fval_all; fval];
                end
                
            else
                for Round = 1:estim.optRound_LPSA %'parfor' will run the parallel computing toolbox
                    k = FalconIC(estim); %initial conditions
                    [xval,fval] = FalconObjFun(estim, k); %objective function
                    %collecting the run's outcome
                    x_all = [x_all; xval];
                    fval_all = [fval_all; fval];
                end
            end
                        
            cost_SA(counter2, counter) = min(fval_all);
            cost_Error_SA(counter2, counter) = std(fval_all);
            params_SA(counter2, counter, find(ismember(Param_original, estim.param_vector))) = ...
                x_all(min(find(ismember(fval_all, min(fval_all)))), :);
            
            if estim.Fast_Option
                if min(fval_all) > CutOff
                    OK = 0;
                end
            end
        end
        
    end
    
    
    % evaluate the identifiability of each parameter
    EvalSA_Temp = cost_SA(:, counter) < CutOff;
    EvalSA = ~EvalSA_Temp;
    Ident_Left = sum(EvalSA(1:floor(length(EvalSA)/2)));
    Ident_Right = sum(EvalSA(ceil(length(EvalSA)/2) + 1:length(EvalSA)));
    if sum([Ident_Left Ident_Right]) == length(EvalSA)-1
        Identifiability = 1;
    elseif Ident_Left || Ident_Right > 0
        Identifiability = 2;
    elseif sum([Ident_Left Ident_Right]) == 0
        Identifiability = 3;
    end
    
    Ident_All(counter)=Identifiability;
    
    try
        % plot the results
        set(0, 'CurrentFigure', thisfig), hold on,
        subplot(NLines, NCols, counter), hold on
        plot(p_SA(:, counter), cost_SA(:, counter), '-.', 'MarkerSize', 5), hold on, 
        
        errorbar(p_SA(:, counter), cost_SA(:, counter), cost_Error_SA(:, counter), 'Color', 'k', 'LineWidth', 1), hold on, 
        plot(p_SA(estim.LPSA_Increments+1, counter), cost_SA(estim.LPSA_Increments+1, counter), 'b*', 'MarkerSize', 5)
        ylabel('MSE')
        title(Param_original(counter))
        hline = refline([0 CutOff]);
        hline.Color = 'r';
        drawnow;
        saveas(thisfig, [estim.FinalFolderName, filesep, 'Identifiability'], 'fig')

    catch
        warning('This LPSA figure is not plotted')
    end
end

try
    suptitle('Local parameter sensitivity analysis');
catch
end

estim = estim_orig;

saveas(thisfig, [estim.FinalFolderName, filesep, 'Identifiability'], 'tif')
saveas(thisfig, [estim.FinalFolderName, filesep, 'Identifiability'], 'fig')
saveas(thisfig, [estim.FinalFolderName, filesep, 'Identifiability'], 'jpg')
saveas(thisfig, [estim.FinalFolderName, filesep, 'Identifiability'], 'svg')
saveas(thisfig, [estim.FinalFolderName, filesep, 'Identifiability'], 'pdf')

%%% Parameter covariance
if size(p_SA, 2) < 20
    disp('generating correlation plot matrix...')
    figure, hold on;
    AllCorrelations = ScatterMatrix(cost_SA, Param_original, ones(size(cost_SA, 1), 1),{}, 1);
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance'], 'tif')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance'], 'fig')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance'], 'jpg')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance'], 'svg')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance'], 'pdf')
    figure, hold on;
    sgtitle('Covariance');
    imagesc(AllCorrelations)
    set(gca, 'xtick', 1:length(Param_original))
    set(gca, 'ytick', 1:length(Param_original))
    set(gca, 'xticklabel', Param_original)
    set(gca, 'yticklabel', Param_original)
    colorbar
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance2'], 'tif')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance2'], 'fig')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance2'], 'jpg')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance2'], 'svg')
    saveas(gca, [estim.FinalFolderName, filesep, 'Covariance2'], 'pdf')
end

estim = estim_orig;

estim.Results.LPSA.ParamNames = estim.param_vector';
estim.Results.LPSA.Identifiability = Ident_All;
estim.Results.LPSA.LPSA_Increments = estim.LPSA_Increments;
estim.Results.LPSA.p_SA = p_SA;
estim.Results.LPSA.cost_SA = cost_SA;
estim.Results.LPSA.CutOff = CutOff;
estim.Results.LPSA.Interpretation = {'1=Identifiable', '2=Partially identifiable', '3=Non-identifiable'};
estim.Results.LPSA.Correlations = AllCorrelations;

Heading = cell(1,2);
Heading(1,1) = {'parameters'};
Heading(1,2) = {'Identifiable'};

LPSA = [estim.param_vector num2cell(estim.Results.LPSA.Identifiability')];

disp('Summary of local parameter sensistivity analysis:')
disp(' ')
disp([Heading; LPSA])
disp(' ')

useexcel = isExcelPresent();
if useexcel        
    xlswrite([estim.FinalFolderName filesep 'Summary_LPSA.xls'], [Heading;LPSA]);
else
    tab = table(estim.param_vector, estim.Results.LPSA.Identifiability', 'VariableNames', Heading);
    writetable(tab, [estim.FinalFolderName, filesep, 'Summary_LPSA.csv'], 'Delimiter', ',');
end

delete('LPSA_TempFile.txt')
end
