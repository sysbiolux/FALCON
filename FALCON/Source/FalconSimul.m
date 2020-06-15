function [estim, StateValueAll, MSE] = FalconSimul(varargin)
% Falcon Simulation
% Runs a Falcon model under specified parameters
% Records the values of states over the course of the simulation
% Arguments: 1) estim (from FalconMakeModel)
% 
% Modified Sebastien De Landtsheer February 2019, sebastien.delandtsheer@uni.lu
% Prof. Thomas Sauter, thomas.sauter@uni.lu

estim = varargin{1};
Plots = 1;
if nargin >1
    Plots = varargin{2};
end
Folder = estim.FinalFolderName; ToSave = 1;

% Initialize necessary parameter and retrieve best parameter set (optParams)
tic
if nargin > 2
    k = varargin{3};
else
    k = estim.Results.Optimization.BestParams;
end

state_names = estim.state_names;
SD = estim.SD;
set(0, 'DefaultTextInterpreter', 'none');

% Perform simulation based on the best parameter set

n = estim.NrStates;
N = numel(estim.Output) - sum(sum(isnan(estim.Output)));

if isfield(estim, 'Weights')
    w = estim.Weights;
else
    w = ones(size(estim.Output));
end

if isfield(estim, 'Lambda')
    l = estim.Lambda;
else
    l = 0;
end

if isfield(estim, 'RegMatrix')
    if isfield(estim.RegMatrix, 'Groups')
        RegGroups = estim.RegMatrix.Groups;
    end
    if isfield(estim.RegMatrix, 'Smooth')
        RegSmooth = estim.RegMatrix.Smooth;
    end
    if isfield(estim.RegMatrix, 'Cluster')
        RegCluster = estim.RegMatrix.Cluster;
    end
end

if isfield(estim, 'Reg')
    if strcmp(estim.Reg, 'none')
        Var = 0;
    elseif strcmp(estim.Reg, 'L1')
        Var = sum(abs(k));
    elseif strcmp(estim.Reg, 'L1Groups')
        Var = sum(sum(abs(k(RegGroups) - repmat(mean(k(RegGroups), 2), 1, size(RegGroups, 2)))));
    elseif strcmp(estim.Reg, 'L1Smooth')
        Var = 0;
        for v = 1:size(RegSmooth, 1)
            Var = Var + sum(abs(k(RegSmooth(v, 2:end)) - k(RegSmooth(v, 1:end-1))));
        end
    elseif strcmp(estim.Reg, 'L2')
        Var = sum(k.^2);
    elseif strcmp(estim.Reg, 'L1/2')
        Var = sum(k.^0.5);
    elseif strcmp(estim.Reg, 'Lx')
        Var = 1/sum(k.^2);
    elseif strcmp(estim.Reg, 'Ldrug')
        Var(1) = sum(k.^0.5);
        Var(2) = 0;
        for v = 1:size(RegGroups, 1)
            km = mean(k(RegGroups(v, :)));
            Var(2) = Var(2) + sum(abs(k(RegGroups(v, :)) - km));
        end
    elseif strcmp(estim.Reg, 'LTriple')
        Var(1) = sum(k.^0.5);
        Var(2) = 0; [S,T] = size(RegCluster);
        kg = sort(k(RegCluster), 2, 'ascend');
        for c = 1:T-1
            for cc = c+1:T
                Var(2) = Var(2) + sum(abs((kg(:, cc) - kg(:, c)) - repmat((cc - c)/T, S, 1)));
            end
        end
        Var(2) = S/Var(2);
        Var(3) = 0;
        for v = 1:size(RegSmooth, 1)
            Var(3) = Var(3) + sum(abs(k(RegSmooth(v, 2:end)) - k(RegSmooth(v, 1:end-1))));
        end
    elseif strcmp(estim.Reg, 'LCluster')
        Var = 0; [S,T] = size(RegCluster);
        kg = sort(k(RegCluster), 2, 'ascend');
        for c = 1:T-1
            for cc = c+1:T
                Var = Var + sum(abs((kg(:, cc) - kg(:, c)) - repmat((cc - c)/T, S, 1)));
            end
        end
        Var = S/Var;
    elseif strcmp(estim.Reg, 'PruneCluster')
        VarA = 0; [S, T] = size(RegCluster);
        kg = sort(k(RegCluster), 2, 'ascend');
        for c = 1:T-1
            for cc = c+1:T
                VarA = VarA + sum(abs((kg(:, cc) - kg(:, c)) - repmat((cc - c)/T, S, 1)));
            end
        end
        Var(2) = S/VarA;
        Var(1) = sum(k.^0.5);            
    end
else
    Var = 0;
end


%initial and successive number of steps for evaluation
%this still needs to be worked on
if n <= 25, initial_t = 10; step_t = 10;
elseif n > 25 && n <= 100, initial_t = 100; step_t = 10;
elseif n > 100, initial_t = 300; step_t = 150;
end

if isfield(estim, 'MaxTime')
    MaxTime = estim.MaxTime;
else
    MaxTime = 0;
end    

ma = estim.ma;
mi = estim.mi;
param_index = estim.param_index;
kmap = k(estim.kInd)'; %extend k including boolean gates. Now k2 has the same length as param_index


pd = param_index(~estim.FixBool, :); %remove fixed Boolean gates from param_index
kA = kmap(find((pd(:, 3) > 0))); %remap k for activating links
kI = kmap(find((pd(:, 4) > 0))); %remap k for inhibiting links
for i = 1:length(estim.IdxInAct) %for each activating link
    ma(estim.IdxOutAct(i), estim.IdxInAct(i)) = kA(i); %map to ma
end
for i = 1:length(estim.IdxInInh) %for each inhibiting link
    mi(estim.IdxOutInh(i), estim.IdxInInh(i)) = kI(i); %map to mi
end

% Get Boolean gate indices in for-loop
BoolMax = max(param_index(:,5)); %number of Boolean gates
if BoolMax > 0 %if there is at least one Boolean gate
    Gate_fill_indices = [];
    Gate_value_fill = [];

    for counter2 = 1:BoolMax %for each Boolean gate: trace the inputs and output
        gate_indices = find(ismember(param_index(:, 5), counter2))';
        current_ma_value = unique(ma(unique(param_index(gate_indices, 1)),param_index(gate_indices, 2)));
        Current_fill_indices = [unique(param_index(gate_indices, 1)), ... %output
            param_index(gate_indices(1), 2), ... %input 1
            param_index(gate_indices(2), 2), ... %input 2
            unique(param_index(gate_indices, 6))]; % type of gate (OR = 1, AND = 2)
        Gate_value_fill = [Gate_value_fill; current_ma_value];
        Gate_fill_indices = [Gate_fill_indices; Current_fill_indices];
    end

    Temp_BoolVal = zeros(n,size(estim.Output_idx, 1));
end

% Load input and output
Input_index = estim.Input_idx;
Output_index = estim.Output_idx;
Inputs = estim.Input;
Measurements = estim.Output;
sds = estim.SD;
sds(sds==0) = 0.05;
sds(isnan(sds))=0.05;
Variances = (sds).^2;

% Evaluation
TimeSoFar = tic;

estim.AllofTheXs = zeros(size(Measurements, 1), 10000, length(state_names));

x = rand(n,size(Measurements, 1)); %initial random values for the nodes
x(Input_index(1, :), :) = Inputs';%(counter_exp,:); %fixing input nodes
cur = 1;
estim.AllofTheXs(:, cur, :) = x';
xmeas = Measurements; %expected values for output nodes
break_point_ss = 1;
runs = initial_t; % Initialize number of time steps

while break_point_ss
    if MaxTime > 0
        if toc(TimeSoFar) > MaxTime
            break
        end
    end
    pre_x = x;
    for counter = 1:runs
        if BoolMax > 0 %calculate values for Boolean gates
            for counter2 = 1:size(Gate_fill_indices, 1)
                if Gate_fill_indices(counter2, 4) == 1 % OR gate

                    Temp_BoolVal(Gate_fill_indices(counter2, 1), :) = ma(Gate_fill_indices(counter2, 1), Gate_fill_indices(counter2, 2)).*...
                        (1-(1-(x(Gate_fill_indices(counter2, 2), :))) .* (1-x(Gate_fill_indices(counter2, 3), :)));

                elseif Gate_fill_indices(counter2,4) == 2 % AND gate

                    Temp_BoolVal(Gate_fill_indices(counter2, 1), :) = ma(Gate_fill_indices(counter2, 1), Gate_fill_indices(counter2, 2)).*...
                        (x(Gate_fill_indices(counter2, 2), :) .* x(Gate_fill_indices(counter2, 3), :));
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%
        %%% core equation %%%
        x = (ma * x) .* (ones(size(x)) - mi * x);
        %%%%%%%%%%%%%%%%%%%%

        if BoolMax > 0 %if there are Boolean gates
            x(Gate_fill_indices(:, 1), :) = Temp_BoolVal(Gate_fill_indices(:, 1), :);
        end
        cur = cur + 1;
        if cur < 10000
            estim.AllofTheXs(:, cur, :) = x';
        end
    end

    if any(sum(abs(pre_x - x)) > estim.SSthresh) %if the network did not reach steady-state
        runs = runs + step_t;
    else %if we are at steady-state
        break_point_ss = 0;
    end

end
xsim = x(Output_index(1, :), :)';
StateValueAll = xsim;
mask = isnan(xmeas);
xsim(mask) = 0; xmeas(mask) = 0;

%calculate the sum-of-squared errors
SSE = (xsim-xmeas).^2 ;
if isfield(estim, 'ObjFunction')
    if strcmp(estim.ObjFunction, 'likelihood')
        MSE = (sum(nansum(SSE ./ Variances .* w))) / N;
    elseif strcmp(estim.ObjFunction, 'MSE')
        MSE = (sum(nansum(SSE .* w))) / N;
    end
else
    MSE = (sum(nansum(SSE .* w)))/N;

end

Diff = MSE + sum(l.*Var);    

Nparams = (sum(k > 0.01));

if isfield(estim, 'Reg')
    if strcmp(estim.Reg, 'L1Groups')
        Std_group = std(k(RegGroups), 0, 2);
        Collapsed = Std_group < 0.01;
        Nparams = sum(Collapsed) + size(RegGroups, 2) * sum(~Collapsed);
    elseif strcmp(estim.Reg, 'L1Smooth')
        Smoothed = (max(k(RegSmooth), [], 2) - min(k(RegSmooth), [], 2)) < 0.1;
        RP =(mean(k(RegSmooth), 2)) > 0.01;
        Nparams = sum(Smoothed .* RP) + size(RegGroups, 2) * sum(~Smoothed .* RP);
    elseif strcmp(estim.Reg, 'Ldrug')
        Std_group = std(k(RegGroups), 0, 2);
        Collapsed = Std_group < 0.01;
        RP=(mean(k(RegGroups), 2)) > 0.01;
        Nparams = sum(Collapsed .* RP) + size(RegGroups, 2) * sum(~Collapsed .* RP);            
    elseif strcmp(estim.Reg, 'LTriple')
            %%% TO DO: check with example

        ParamIdentity = (abs(bsxfun(@minus, k, k')) < 0.01);
        ParamMatters = zeros(length(k));
        for kline = 1:size(RegCluster,1)
            for kcol1 = 1:size(RegCluster,2) - 1
                for kcol2 = kcol1 + 1:size(RegCluster, 2)
                    ParamMatters(RegCluster(kline, kcol1), RegCluster(kline, kcol2)) = 1;
                end
            end
        end
        for kline = 1:size(RegSmooth, 1)
            for kcol1 = 1:size(RegSmooth, 2) - 1
                for kcol2 = kcol1 + 1:size(RegSmooth, 2)
                    ParamMatters(RegSmooth(kline, kcol1), RegSmooth(kline, kcol2)) = 1;
                end
            end
        end
        ParamMatters = ParamMatters .* (~eye(length(k)));
        ParamSignif = k' > 0.01;
        NCollapsed = ceil(sum((sum(ParamMatters.*ParamIdentity, 2) > 0) .* ParamSignif));
        Nparams = estim.NrParams - NCollapsed;
            %%% TO DO: check with example

    elseif strcmp(estim.Reg, 'LCluster')
        Par = (k(RegCluster));
        Par = sort(Par, 2);
        Dist = Par(:, 2:end) - Par(:, 1:end-1);
        Collapsed = Dist < 0.01;
        Nparams = size(RegCluster, 1) + sum(~Collapsed(:));
    elseif strcmp(estim.Reg, 'PruneCluster')
        Par = (k(RegCluster));
        Par = sort(Par, 2);
        Sig = Par > 0.01;
        Dist = Par(:, 2:end) - Par(:, 1:end-1);
        Collapsed = Dist < 0.01;
        Sig = Sig .* [ones(size(Collapsed, 1), 1), ~Collapsed];
        Nparams = sum(Sig(:));
    end
end

AIC = N .* log(MSE) + 2 .* Nparams;
BIC = N .* log(MSE) + Nparams * log(N);

fprintf('MSE= %d \t reg cost= %d \t total= %d \t BIC= %d \n', MSE, sum(l.*Var), Diff, BIC);

if Plots
    if estim.PlotFitSummary
        % Plot simulated state value mapped on experimental data
        num_plots = size(estim.Output, 2);
        NLines = ceil(sqrt(num_plots));
        NCols = ceil(num_plots/NLines);

        % Plot molecular profiles
        h1 = figure; hold on
        for counter = 1:num_plots
            subplot(NLines, NCols, counter), hold on,

            % Plot experimental data first (in green)
            if ~isempty(SD)
                errorbar(1:size(Measurements, 1),Measurements(:, counter),SD(:, counter), 'gs', 'LineWidth', 1, 'MarkerSize', 2, 'Color', [0.4 0.6 0]), hold on,
            else
                errorbar(1:size(Measurements, 1), Measurements(:, counter), zeros(size(Measurements, 1), 1), 'gs', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [0.4 0.6 0]), hold on,
            end

            % Plot simulated data on top
            plot(1:size(Measurements, 1),StateValueAll(:, counter), 'b.', 'MarkerSize', 20 / sqrt(num_plots))

            % Figure adjustment
            axis([0 size(Measurements,1)+1 0 1.1])
            set(gca, 'XTick', 1:length(estim.Annotation))
            set(gca, 'XTickLabel', estim.Annotation)
            set(gca, 'XTickLabelRotation', 45)
            set(gca, 'fontsize', 15/sqrt(num_plots))
            t = title(state_names(Output_index(1, counter)));
            xt = xlabel('experimental condition');
            set(xt, 'fontsize', 15/sqrt(num_plots))
            set(t, 'fontsize', 25/sqrt(num_plots))
            hold off
        end
        if ToSave
            saveas(h1, [Folder,'\Fitting_plot'], 'tif')
            saveas(h1, [Folder,'\Fitting_plot'], 'pdf')
        end
    end

    if estim.PlotFitIndividual
        % Plot simulated state value mapped on experimental data
        % Plot molecular profiles
        for counter = 1:length(state_names)
            h1 = figure; hold on

            % Plot experimental data first (in green)
            for counter_plot = 1:size(Output_index(1, :), 2)
                current_counter_plot = Output_index(1, :);
                if counter == current_counter_plot(counter_plot) %Panuwat's logic ;-)
                    if ~isempty(SD)
                        errorbar(1:size(Measurements, 1), Measurements(:, counter_plot), SD(:, counter_plot), 'gs', 'LineWidth', 1, 'MarkerSize', 2, 'Color', [0.4 0.6 0]), hold on,
                    else
                        errorbar(1:size(Measurements, 1), Measurements(:, counter_plot), zeros(size(Measurements, 1), 1), 'gs', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [0.4 0.6 0]), hold on,
                    end
                end
            end

            % Plot simulated data on top
            plot(1:size(Measurements,1), StateValueAll(counter, :), 'b.', 'MarkerSize', 15)

            % Figure adjustment
            axis([0 size(Measurements,1)+1 0 1.1])
            set(gca, 'XTick', 1:length(estim.Annotation))
            set(gca, 'XTickLabel', estim.Annotation)
            set(gca, 'XTickLabelRotation', 45)
            set(gca, 'fontsize', 15)
            t = title(state_names(counter));
            xt = xlabel('exp');
            y = ylabel('state-value');
            set(xt,'fontsize', 15)
            set(y, 'fontsize', 15)
            set(t, 'fontsize', 15)
            hold off
            if ToSave
                saveas(h1, [Folder, '\', cell2mat(state_names(counter)) '_plot'], 'tif')
                saveas(h1, [Folder, '\', cell2mat(state_names(counter)) '_plot'], 'pdf')
            end

        end
    end

    if estim.PlotHeatmapCost && sum(std(estim.Output_idx)) == 0
        % Plot optimal cost for each experiment
        MapDiffs = [SSE; mean(SSE, 1)]; MapDiffs = [MapDiffs, mean(MapDiffs, 2)];
    %     hm = HeatMap(MapDiffs, 'RowLabels', cellstr([estim.Annotation; 'mean']), 'ColumnLabels', [estim.state_names(estim.Output_idx(1, :)), 'mean'], 'Colorbar', 'on', 'Colormap', hot, 'Symmetric', false);
    %     addTitle(hm, 'Cross-error Analysis: Heatmap');
        figure, hold on
        imagesc(MapDiffs), colorbar,
        set(gca, 'xtick', 1:size(MapDiffs, 2))
        set(gca, 'xticklabel', [estim.state_names(estim.Output_idx(1, :)), 'mean'])
        set(gca, 'xticklabelrotation', 90)
        set(gca, 'ytick', 1:size(MapDiffs, 1))
        set(gca, 'yticklabel', cellstr([estim.Annotation; 'mean']))
        axis([0.5, size(MapDiffs, 2)+0.5, 0.5, size(MapDiffs, 1)+0.5])
        colormap('hot')
        title('Cross-error Analysis: Heatmap')
        if ToSave
    %         fighm = plot(hm);
            saveas(gcf,[Folder, '\CrossErrorHeatMap'], 'tif');
            saveas(gcf,[Folder, '\CrossErrorHeatMap'], 'pdf');
        end
        ploty = x'; ploty = ploty(:, estim.Output_idx(1, :)); ploty = ploty(:);
        plotx = Measurements(:);

        figure, cp = plot(plotx,ploty, '.k'); xlabel('Measured'), ylabel('Simulated');
        mdl = fitlm(plotx,ploty); hold on
        R2 = mdl.Rsquared.Adjusted; R2 = (round(R2*1000))/1000;
        Pval = mdl.Coefficients{2,4};
        a = mdl.Coefficients{1, 1}; b = mdl.Coefficients{2, 1};
        plot([0 1], [a, a+b], '--k')
        title(['Simulated versus Measurements: R^2 = ', num2str(R2),', p = ', num2str(Pval) ]);
        set(gca, 'TickLabelInterpreter', 'latex')
        axis([-0.05 1.05 -0.05 1.05])
        if ToSave
            saveas(cp, [Folder, '\CorrelationPlot'], 'tif');
            saveas(cp, [Folder, '\CorrelationPlot'], 'pdf');
        end
    end

    if estim.PlotStateSummary
        h4 = figure; hold on
        NrExps = length(estim.Output(:, 1));
        for i = 1:NrExps %for each exp
            subplot(NrExps, 1, i)
            Y = state_names; X = Y;
            set(gca, 'xtick', 1:length(X))
            set(gca, 'xticklabel', X)
            hold on
            Y2mean = x(:, i);
            Y2std = zeros(size(Y2mean));    
            errorbar(Y2mean, Y2std, '.b', 'LineWidth', 1)
            hold on
            if ~isempty(SD)
                errorbar(Output_index(i, :), Measurements(i, :), SD(i, :), '.g', 'Color', [0.4 0.6 0])
            else
                errorbar(Output_index(i, :), Measurements(i, :), zeros(size(Measurements(i,:))), '.g', 'Color', [0.4 0.6 0])
            end
            axis([0, length(X)+1, 0, 1.1])
        end
        hold off
        if ToSave
            saveas(h4, [Folder, '\Measured_vs_Simul'], 'tif')
            saveas(h4, [Folder, '\Measured_vs_Simul'], 'pdf')
        end
    end

    if estim.PlotStateEvolution
        T = 0;
        for exp = 1:size(Measurements, 1)
            for sim = 2:size(estim.AllofTheXs, 2)
                if (sum(abs(estim.AllofTheXs(exp, sim, :) - estim.AllofTheXs(exp, sim-1, :)))) < 0.01
                    T = max(sim, T);
                    break
                end
            end
        end


        num_plots = length(Output_index);
        NLines = ceil(sqrt(num_plots));
        NCols = ceil(num_plots/NLines);

        h51 = figure; hold on,
         suptitle('Dynamics through the simulation (outputs)');
        for p = 1:size(Output_index, 1)
            subplot(NLines, NCols, p)
            plot(squeeze(estim.AllofTheXs(p, 1:T, Output_index(p, :))))
            axis([0, T+1, 0, 1]);
            hold off
        end
        legend(state_names(Output_index(p, :)))

        if ToSave
            saveas(h51, [Folder, '\ConvergenceOutputNodes'], 'tif')
            saveas(h51, [Folder, '\ConvergenceOutputNodes'], 'pdf')
        end

        h52 = figure; hold on,

        suptitle('Dynamics through the simulation (all nodes)')
        for p = 1:size(Output_index, 1)
            subplot(NLines, NCols, p)
            plot(squeeze(estim.AllofTheXs(p, 1:T, :)))
            axis([0, T+1, 0, 1]);
            hold off
        end
          legend(state_names(:))
        if ToSave
            saveas(h52, [Folder, '\ConvergenceAllNodes'], 'tif')
            saveas(h52, [Folder, '\ConvergenceAllNodes'], 'pdf')
        end
    end
end
StateValueAll = x';
estim.Results.Optimization.StateValueAll = StateValueAll;
estim.Results.Optimization.SquaredErrors = SSE;
end

