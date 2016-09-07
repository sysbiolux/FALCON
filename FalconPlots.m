function FalconPlots(estim)
% FalconPlots displays all the different graphs at the end of the optimization
% FalconPlots(estim)
%
% :: Input values ::
% estim        complete model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

MeanStateValueAll = estim.MeanStateValueAll;
Output_index=estim.Output_idx;
state_names=estim.state_names;
Measurements=estim.Output;
SD=estim.SD;
bestx = estim.bestx;

%for states
num_plots=size(estim.Output,2);
NLines=ceil(sqrt(num_plots));
NCols=ceil(num_plots/NLines);
%for parameters
num_plots_params=size(estim.param_vector,1);
NLines_params=ceil(sqrt(num_plots_params));
NCols_params=ceil(num_plots_params/NLines_params);

%% Plot Fit Evolution
if isfield(estim.Results,'FitEvol')
    PlotCosts =  estim.Results.FitEvol.PlotCosts;
    Cost1 = estim.Results.FitEvol.Cost1;
    Cost2 = estim.Results.FitEvol.Cost2;
    Cost3 = estim.Results.FitEvol.Cost3;
    
    hconv=figure;
    plot(1:max([length(Cost1),length(Cost2),length(Cost3)]),PlotCosts)
    title('SSE convergence during optimization')
    xlabel('iteration'); ylabel('SSE');
end

%% Plot simulated state value mapped on experimental data
if isfield(estim.Results,'Optimisation')
    Diffs= estim.Results.Optimisation.Diffs;
    StdStateValueAll = estim.Results.Optimisation.StdStateValueAll;
    FittingCost = estim.Results.Optimisation.FittingCost;
    
    % show plotting evolution (needed to show that the model really ended up in a global optimum)
    h11=  figure;
    plot(sort(FittingCost), 'o-');
    xlabel('number of fits');
    ylabel('goodnes s of fits');
    title('Fitting cost');
    
    
    % if graphs_opt(1)
    h1=figure; hold on
    for counter=1:num_plots
        subplot(NLines,NCols,counter), hold on,
        
        % Plot experimental data first (in green)
        
        if ~isempty(SD)
            errorbar(1:size(Measurements,1),Measurements(:,counter),SD(:,counter),'gs','LineWidth',3,'MarkerSize',5), hold on,
        else
            errorbar(1:size(Measurements,1),Measurements(:,counter),zeros(size(Measurements,1),1),'gs','LineWidth',3,'MarkerSize',5), hold on,
        end
        
        % Plot simulated data on top (error bar in red, mean in blue)
        errorbar(1:size(Measurements,1),MeanStateValueAll(:,Output_index(1,counter)),StdStateValueAll(:,Output_index(1,counter)),'r.','LineWidth',3), hold on,
        plot(1:size(Measurements,1),MeanStateValueAll(:,Output_index(1,counter)),'b*','MarkerSize',25/sqrt(num_plots))
        
        % Figure adjustment
        axis([0 size(Measurements,1)+1 0 1.21])
        set(gca,'fontsize',25/sqrt(num_plots))
        t=title(state_names(Output_index(1,counter)));
        x=xlabel('exp');
        y=ylabel('state-value');
        set(x,'fontsize',25/sqrt(num_plots))
        set(y,'fontsize',25/sqrt(num_plots))
        set(t,'fontsize',35/sqrt(num_plots))
        hold off
    end
    % end
    
    % Plot simulated state value mapped on experimental data
    % Plot molecular profiles
    % if graphs_opt(2)
    for counter=1:length(state_names)
        h1=figure; hold on
        
        % Plot experimental data first (in green)
        for counter_plot=1:size(Output_index(1,:),2)
            current_counter_plot=Output_index(1,:);
            if counter==current_counter_plot(counter_plot)
                if ~isempty(SD)
                    errorbar(1:size(Measurements,1),Measurements(:,counter_plot),SD(:,counter_plot),'gs','LineWidth',3,'MarkerSize',5), hold on,
                else
                    errorbar(1:size(Measurements,1),Measurements(:,counter_plot),zeros(size(Measurements,1),1),'gs','LineWidth',3,'MarkerSize',5), hold on,
                end
            end
        end
        
        % Plot simulated data on top (error bar in red, mean in blue)
        errorbar(1:size(Measurements,1),MeanStateValueAll(:,counter),StdStateValueAll(:,counter),'r.','LineWidth',3), hold on,
        plot(1:size(Measurements,1),MeanStateValueAll(:,counter),'b*','MarkerSize',15)
        
        
        % Figure adjustment
        axis([0 size(Measurements,1)+1 0 1.21])
        set(gca,'fontsize',15)
        t=title(state_names(counter));
        set(t,'interpreter','none')
        x=xlabel('exp');
        y=ylabel('state-value');
        set(x,'fontsize',15)
        set(y,'fontsize',15)
        set(t,'fontsize',15)
        hold off
        
    end
    % end
    
    % PlotHeatmapCost: heatmaps of optimal costs for each output for each experimental condition
    % if graphs_opt (3)
    figure;
    hm=imagesc(Diffs);
    colorbar
    X_axis_name=estim.state_names(estim.Output_idx(1,:));
    
    set(gca,'XTick',1:size(estim.Output_idx,1))
    set(gca,'xticklabel',X_axis_name)
    colormap('hot')
    title('Cross-error Analysis: Heatmap')
    ylabel('Experiments')
    % end
    
    % PlotStateSummary: graph of state values at steady-state (all-in-1)
    % if graphs_opt(5)
    h4=figure; hold on
    NrExps=length(estim.Output(:,1));
    if NrExps > 5
        NrExps = 5;
    end
    for i=1:NrExps %for each exp
        subplot(NrExps,1,i)
        Y=state_names; X=Y;
        set(gca, 'xtick', 1:length(X))
        set(gca,'xticklabel',X)
        hold on
        if ~isempty(SD)
            errorbar(Output_index(i,:), Measurements(i,:), SD(i,:),'gs','LineWidth',3,'MarkerSize',5)
        else
            errorbar(Output_index(i,:), Measurements(i,:), zeros(size(Measurements(i,:))),'gs','LineWidth',3,'MarkerSize',5)
        end
        hold on
        Y2mean=MeanStateValueAll(i,:);
        Y2std=StdStateValueAll(i,:);
        errorbar(Y2mean, Y2std,'r.','LineWidth',3)
        plot(1:length(state_names),Y2mean,'b*','MarkerSize',10)
        axis([0, length(X)+1, 0, 1])
        ylabel(['Exp' num2str(i)])
    end
    hold off
    if NrExps > 5
        suptitle('Compared Simulated values and Measurement (all states)')
    else
        suptitle('Examples of compared Simulated values and Measurement (all states)')
    end
    
    % PlotStateEvolution: graph of state values over the course of the simulation (two graphs: all-in-1 and individual)
    % if graphs_opt (6)
    T=0;
    for exp=1:size(Measurements,1)
        for sim=2:size(estim.AllofTheXs,2)
            if (sum(abs(estim.AllofTheXs(exp,sim,:)-estim.AllofTheXs(exp,sim-1,:))))<0.01
                T=max(sim,T);
                break
            end
        end
    end
    
    h51=figure; hold on,
    
    NrExps=length(estim.Output(:,1));
    if NrExps > 5
        NrExps = 5;
    end
    for p=1:NrExps %for each exp
        subplot(NrExps,1,p)
        plot(squeeze(estim.AllofTheXs(p,1:T,Output_index(p,:))))
        legend(state_names(Output_index(p,:)),'Location','EastOutside')
        ylabel(['Exp' num2str(p)])
        axis([0, T+1, 0, 1]);
        hold off
    end
    if NrExps > 5
        suptitle('Dynamics through the simulation (outputs)')
    else
        suptitle('Examples of dynamics through the simulation (outputs)')
    end
    
    h52=figure; hold on,
    NrExps=length(estim.Output(:,1));
    if NrExps > 5
        NrExps = 5;
    end
    for p=1:NrExps %for each exp
        subplot(NrExps,1,p)
        plot(squeeze(estim.AllofTheXs(p,1:T,:)))
        legend(state_names(:),'Location','EastOutside')
        ylabel(['Exp' num2str(p)])
        axis([0, T+1, 0, 1]);
        hold off
    end
    if NrExps > 5
        suptitle('Dynamics through the simulation (all nodes)')
    else
        suptitle('Examples of dynamics through the simulation (all nodes)')
    end
    
end

%% Resampling
if isfield(estim.Results,'Resampling')
    
    Costs = estim.Results.Resampling.Costs;
    Ks = estim.Results.Resampling.OptimisedParameter;
    % if graphs_Resampling(1)
    disp(['Cost for perturbated measurements: ', num2str(mean(Costs)), ' +- ', num2str(std(Costs))])
    h1= figure; hold on;
    errorbar(mean(Ks),std(Ks),'.b','LineWidth',2); hold on;
    ylim([0 1])
    plot(bestx,'*r','MarkerSize',12);
    title('Distribution of optimised parameters after resampling')
    hold on;
    set(gca,'Xtick',1:length(Ks),'XTickLabel',estim.param_vector,'XGrid','on');
    legend('RangeResampling','BestParam')
end

%% LPSA
if isfield(estim.Results,'LPSA')
    
    CutOff = estim.Results.LPSA.CutOff;
    p_SA= estim.Results.LPSA.p_SA;
    Param_original=estim.param_vector;
    cost_SA= estim.Results.LPSA.cost_SA;
    p_increment= estim.Results.LPSA.p_increment;
    % if graphs_LPSA (1)
    figure,
    for counter = 1:length(estim.param_vector)
        subplot(NLines_params, NCols_params, counter),
        plot(p_SA(:,counter),cost_SA(:,counter))
        hold on, plot(p_SA(p_increment+1,counter),cost_SA(p_increment+1,counter),'b*','MarkerSize',15)
        min_val=min(cost_SA(:,counter));
        min_index=find(min_val==cost_SA(:,counter));
        hold on, plot(p_SA(min_index,counter),cost_SA(min_index,counter),'rs','MarkerSize',15)
        
        xlabel('Parameter range')
        ylabel('SSE')
        title(Param_original(counter))
        hline=refline([0 CutOff]);
        hline.Color = 'r';
    end
    suptitle('Local parameter sensitivity analysis');
end
%% KO
if isfield(estim.Results,'KnockOut')
    Xtitles = estim.Results.KnockOut.Parameters;
    AIC_merge = estim.Results.KnockOut.AIC_values;
    
    figure; hold on;
    for counter = 1: length(AIC_merge);
        h=bar(counter, AIC_merge(counter));
        if AIC_merge(counter) < AIC_merge(1);
            set(h,'FaceColor','g');
        elseif AIC_merge(counter) == AIC_merge(1);
            set(h,'FaceColor','b');            
        else AIC_merge(counter);
            set(h,'FaceColor','k');
        end
    end
    hold off
    hline=refline([0 AIC_merge(1)]);
    hline.Color = 'r';
    
    set(gca,'XTick',1:length(AIC_merge)+1)
    set(gca, 'XTicklabel', Xtitles);
    title('AIC');
    xlabel('');
    ylabel('Akaike Information Criterion (AIC)');
    hold off
    
end
end