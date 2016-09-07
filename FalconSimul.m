function [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(varargin)
% FalconSimul resimulates the optimized network and displays the summarised results as plots
% [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

% :: Input values ::
% estim                     complete model definition
% bestx                     the vector of best optimised parameters
% 5-values logical vector   0= do not show plot, 1= show plot
%   - PlotFitSummary        graph of state values at steady-state versus measurements (all-in-1)
%   - PlotFitIndividual     graph of state values at steady-state versus measurements (individual)
%   - PlotHeatmapCost       heatmaps of optimal costs for each output for each experimental condition
%   - PlotStateSummary      graph of state values at steady-state (all-in-1)
%   - PlotStateEvolution    graph of state values over the course of the simulation (two graphs: all-in-1 and individual)
% FinalFolderName           name of the folder for saving results
%
% :: Output values :: 
% MeanStateValueAll         mean states values
% StdStateValueAll          standard deviation of state values
% MeanCostAll               mean fitting cost (SSE)
% StdCostAll                standard deviation of fitting cost
% estim                     updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

estim=varargin{1}; optParams=varargin{2}; graphs=varargin{3};
ToSave=0;
if nargin>3, Folder=varargin{4}; ToSave=1; end

% Initialize necessary parameter and retrieve best parameter set (optParams)
tic
k=optParams;
MeanStateValueAll=[];
StdStateValueAll=[];
MeanCostAll=[];
StdCostAll=[];
Diffs=[];
Max=10000;
Input_index=estim.Input_idx;
Output_index=estim.Output_idx;
Inputs=estim.Input;
state_names=estim.state_names;
Measurements=estim.Output;
SD=estim.SD;
param_index=estim.param_index;
ma=estim.ma;
mi=estim.mi;

% Perform simulation based on the best parameter set

n=estim.NrStates;

%initial and successive number of steps for evaluation
if n<=25, initial_t=10; step_t=10;
elseif n>25 && n<=100, initial_t=100; step_t=10;
elseif n>100, initial_t=300; step_t=150;
end

k2=k(estim.kInd)'; %extend k including boolean gates. Now k2 has the same length as param_index

pd=param_index(~estim.FixBool,:); %shorter name for param_index
kA=k2(find((pd(:,3)>0))); %remap k for activating links
kI=k2(find((pd(:,4)>0))); %remap k for inhibiting links
for i=1:length(estim.IdxInAct) %for each activating link
    ma(estim.IdxOutAct(i),estim.IdxInAct(i))=kA(i); %map to ma
end
for i=1:length(estim.IdxInInh) %for each inhibiting link
    mi(estim.IdxOutInh(i),estim.IdxInInh(i))=kI(i); %map to mi
end
%%%


% Get Boolean gate indices in for-loop
BoolMax=max(param_index(:,5));
if BoolMax>0
    Gate_fill_indices=[];
    Gate_value_fill=[];

    for counter2=1:BoolMax
        gate_indices=find(ismember(param_index(:,5),counter2))';
        current_ma_value=unique(ma(unique(param_index(gate_indices,1)),param_index(gate_indices,2)));
        Current_fill_indices=[unique(param_index(gate_indices,1)), ... %output
            param_index(gate_indices(1),2), ... %input 1
            param_index(gate_indices(2),2), ... %input 2
            unique(param_index(gate_indices,6))]; % type of gate (OR = 1, AND = 2)
        Gate_value_fill=[Gate_value_fill; current_ma_value];
        Gate_fill_indices=[Gate_fill_indices; Current_fill_indices];
    end

    Temp_BoolVal=zeros(1,BoolMax);
end

% Note: it is possible to perform multiple simulations and collect results here
estim.AllofTheXs=zeros(size(Measurements,1),Max,length(state_names));

for counter_exp=1:size(Measurements,1)
    Collect_x=[];
    diff_ALL=[];
    x=rand(1,size(ma,2))';
    x(Input_index(counter_exp,:))=Inputs(counter_exp,:);
    cur=1;
    estim.AllofTheXs(counter_exp,cur,:)=x';
    xmeas=Measurements(counter_exp,:);
        
    idx_NaN=find(isnan(xmeas));
        
    diff=0; % Initialize fitting cost
    break_point_ss=1;

    runs=initial_t; % Initialize number of time steps

    while  break_point_ss
        pre_x=x;
        for counter=1:runs
            if BoolMax>0
                for counter2=1:size(Gate_fill_indices,1)
                    if Gate_fill_indices(counter2,4)==1 % OR gate

                        Temp_BoolVal(counter2)=ma(Gate_fill_indices(counter2,1),Gate_fill_indices(counter2,2))*...
                            (1-(1-(x(Gate_fill_indices(counter2,2))))*(1-x(Gate_fill_indices(counter2,3))));

                    elseif Gate_fill_indices(counter2,4)==2 % AND gate

                        Temp_BoolVal(counter2)=ma(Gate_fill_indices(counter2,1),Gate_fill_indices(counter2,2))*...
                            (x(Gate_fill_indices(counter2,2))*x(Gate_fill_indices(counter2,3)));
                    end
                end
            end

            x=(ma*x).*(ones(numel(x),1)-mi*x);

            if BoolMax>0
                for counter3=1:size(Gate_fill_indices,1)
                    x(Gate_fill_indices(counter3,1))=Temp_BoolVal(counter3);
                end
            end                                
            cur=cur+1;
            if cur<Max
                estim.AllofTheXs(counter_exp,cur,:)=x';
            end
        end
        if sum(abs(pre_x-x))>estim.SSthresh
            runs=runs+step_t;
        else
            break_point_ss=0;
        end

    end

    if ~isempty(idx_NaN)
        OutNaN=x(Output_index(counter_exp,idx_NaN));
        x(Output_index(counter_exp,idx_NaN))=0;
        xmeas(idx_NaN)=0;
    end

    diff=diff+sum((x(Output_index(counter_exp,:))'-xmeas).^2);
    Collect_x=[Collect_x; x'];
    diff_ALL=[diff_ALL; diff];
    
    % Collect state values in each round of simulation
    
    if ~isempty(idx_NaN)
        for counter=1:length(idx_NaN)
            Collect_x(:,Output_index(counter_exp,idx_NaN(counter)))=OutNaN(counter);
        end
    end
    
    MeanAllState=mean(Collect_x,1);
    StdAllState=std(Collect_x,0,1);
    Diffs=[Diffs;(x(Output_index(counter_exp,:))'-xmeas).^2];
    mean_diff_ALL=mean(diff_ALL,1);
    std_diff_ALL=std(diff_ALL,0,1);
    
    MeanStateValueAll=[MeanStateValueAll; MeanAllState];
    StdStateValueAll=[StdStateValueAll; StdAllState];
    MeanCostAll=[MeanCostAll; mean_diff_ALL];
    StdCostAll=[StdCostAll; std_diff_ALL];
    
end

estim.Results.Optimisation.BestStates = Collect_x;

% PlotFitSummary: graph of state values at steady-state versus measurements (all-in-1)

if graphs(1)
    % Plot simulated state value mapped on experimental data
    num_plots=size(estim.Output,2);
    NLines=ceil(sqrt(num_plots));
    NCols=ceil(num_plots/NLines);

    % Plot molecular profiles
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
    if ToSave
        saveas(h1,[Folder,filesep,'Fitting_plot'],'tif')
        saveas(h1,[Folder,filesep,'Fitting_plot'],'fig')

    end
end

% PlotFitIndividual: graph of state values at steady-state versus measurements (individual)

if graphs(2)
    % Plot simulated state value mapped on experimental data
    % Plot molecular profiles
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
        if ToSave
            saveas(h1,[Folder, filesep, 'State_' cell2mat(state_names(counter)) '_plot'],'tif')
            saveas(h1,[Folder, filesep, 'State_' cell2mat(state_names(counter)) '_plot'],'fig')

        end

    end
end

% PlotHeatmapCost: heatmaps of optimal costs for each output for each experimental condition

if graphs(3) && sum(std(estim.Output_idx))==0    
    figure;
    hm=imagesc(Diffs);
    colorbar
    X_axis_name=estim.state_names(estim.Output_idx(1,:));
    
    set(gca,'XTick',1:size(estim.Output_idx,1))
    set(gca,'xticklabel',X_axis_name)
    colormap('hot')
    title('Cross-error Analysis: Heatmap')
    ylabel('Experiments')
    
    if ToSave
        fighm=hm;
        saveas(fighm,[Folder, filesep, 'CrossErrorHeatMap'],'tif');
        saveas(fighm,[Folder, filesep, 'CrossErrorHeatMap'],'fig');
    end
%     figure, hist(Diffs(:)); title('Cross-error Analysis: Histogram')
%     if ToSave
%         saveas(gcf,[Folder, filesep, 'CrossErrorHistogram'],'tif')
%         saveas(gcf,[Folder, filesep, 'CrossErrorHistogram'],'fig')
%     end
end

% PlotStateSummary: graph of state values at steady-state (all-in-1)

if graphs(4)
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

    
    if ToSave
        saveas(h4,[Folder, filesep, 'Measured_vs_Simul'],'tif')
        saveas(h4,[Folder, filesep, 'Measured_vs_Simul'],'fig')

    end
end

% PlotStateEvolution: graph of state values over the course of the simulation (two graphs: all-in-1 and individual)

if graphs(5)
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
    
    if ToSave
        saveas(h51,[Folder, filesep, 'ConvergenceOutputNodes'],'tif')
        saveas(h51,[Folder, filesep, 'ConvergenceOutputNodes'],'fig')

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
    
    if ToSave
        saveas(h52,[Folder, filesep,'ConvergenceAllNodes'],'tif')
        saveas(h52,[Folder, filesep,'ConvergenceAllNodes'],'fig')

    end
end

estim.Results.Optimisation.Diffs = Diffs
estim.Results.Optimisation.StdStateValueAll = StdStateValueAll
end