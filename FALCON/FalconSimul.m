function [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(varargin)
% Falcon Simulation
% Runs a Falcon model under specified parameters
% Records the values of states over the course of the simulation
% Arguments: 1) estim (from FalconMakeModel)
% 2) the vector of optimal parameters
% 3) 5-values logical vector:
% graph of state values at ss versus measurements (all in 1)
% graph of state values at ss versus measurements (individual)
% heatmaps of optimal costs for each output for each condition absolute cost
% graph of state values at ss (all in 1)
% graph of state values over the course of the simulation (two graphs)
% Modified Sebastien De Landtsheer July 2016,
% sebastien.delandtsheer@uni.lu
% Panuwat Trairatphisan, University of Luxembourg, 10/2014, panuwat.trairatphisan@uni.lu

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

num_plots=length(Output_index);
NLines=ceil(sqrt(num_plots));
NCols=ceil(num_plots/NLines);
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

    Temp_BoolVal=zeros(n,size(estim.Output_idx,1));
end

% Note: it is possible to perform multiple simulations and collect results here
estim.AllofTheXs=zeros(size(Measurements,1),Max,length(state_names));

% Collect_x=[];
% diff_ALL=[];
x=rand(n,size(Measurements,1)); %initial random values for the nodes
x(Input_index(1,:),:)=Inputs';
cur=1;
estim.AllofTheXs(:,cur,:)=x';
xmeas=Measurements;

idx_NaN=find(isnan(xmeas));

diff=0; % Initialize fitting cost
break_point_ss=1;

runs=initial_t; % Initialize number of time steps

while  break_point_ss
    pre_x=x;
    for counter=1:runs
        if BoolMax>0 %calculate values for Boolean gates
            for counter2=1:size(Gate_fill_indices,1)
                if Gate_fill_indices(counter2,4)==1 % OR gate

                    Temp_BoolVal(Gate_fill_indices(counter2,1),:)=ma(Gate_fill_indices(counter2,1),Gate_fill_indices(counter2,2)).*...
                        (1-(1-(x(Gate_fill_indices(counter2,2),:))).*(1-x(Gate_fill_indices(counter2,3),:)));

                elseif Gate_fill_indices(counter2,4)==2 % AND gate

                    Temp_BoolVal(Gate_fill_indices(counter2,1),:)=ma(Gate_fill_indices(counter2,1),Gate_fill_indices(counter2,2)).*...
                        (x(Gate_fill_indices(counter2,2),:).*x(Gate_fill_indices(counter2,3),:));
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%
        %%% core equation %%%
        x=(ma*x).*(ones(size(x))-mi*x);
        %%%%%%%%%%%%%%%%%%%%

        if BoolMax>0 %if there are Boolean gates
            x(Gate_fill_indices(:,1),:)=Temp_BoolVal(Gate_fill_indices(:,1),:);
        end                                
        cur=cur+1;
        if cur<Max
            estim.AllofTheXs(:,cur,:)=x';
        end
    end
    if any(sum(abs(pre_x-x))>estim.SSthresh) %if the network did not reach steady-state
        runs=runs+step_t;
    else %if we are at steady-state
        break_point_ss=0;
    end

end

xsim=x(Output_index(1,:),:)';
mask=isnan(xmeas);
xsim(mask)=0; xmeas(mask)=0;

%calculate the sum-of-squared errors
diff=sum(sum((xsim-xmeas).^2));

disp(diff)

diff_ALL=diff;

MeanAllState=xsim;
StdAllState=zeros(size(xsim));
Diffs=(xsim-xmeas).^2;

MeanStateValueAll=[MeanStateValueAll; MeanAllState];
StdStateValueAll=[StdStateValueAll; StdAllState];
MeanCostAll=diff;
StdCostAll=0;


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
            errorbar(1:size(Measurements,1),Measurements(:,counter),SD(:,counter),'gs','LineWidth',1,'MarkerSize',2,'Color',[0.4 0.6 0]), hold on,
        else
            errorbar(1:size(Measurements,1),Measurements(:,counter),zeros(size(Measurements,1),1),'gs','LineWidth',1,'MarkerSize',1,'Color',[0.4 0.6 0]), hold on,
        end

        % Plot simulated data on top
        plot(1:size(Measurements,1),MeanStateValueAll(:,counter),'b.','MarkerSize',20/sqrt(num_plots))

        % Figure adjustment
        axis([0 size(Measurements,1)+1 0 1.1])
        set(gca,'XTick', [1:length(estim.Annotation)])
        set(gca,'XTickLabel', estim.Annotation)
        xtickangle(45)
        set(gca,'fontsize',15/sqrt(num_plots))
%         set(gca,'XMinorGrid','on')
        t=title(state_names(Output_index(1,counter)));
        xt=xlabel('experimental condition');
        yt=ylabel('state-value');
        set(xt,'fontsize',15/sqrt(num_plots))
        set(t,'fontsize',25/sqrt(num_plots))
        hold off
    end
    if ToSave
        saveas(h1,[Folder,'\Fitting_plot'],'tif')
    end
end

if graphs(2)
    % Plot simulated state value mapped on experimental data
    % Plot molecular profiles
    for counter=1:length(state_names)
        h1=figure; hold on

        % Plot experimental data first (in green)
        for counter_plot=1:size(Output_index(1,:),2)
            current_counter_plot=Output_index(1,:);
            if counter==current_counter_plot(counter_plot) %Panuwat's logic
                if ~isempty(SD)
                    errorbar(1:size(Measurements,1),Measurements(:,counter_plot),SD(:,counter_plot),'gs','LineWidth',1,'MarkerSize',2, 'Color',[0.4 0.6 0]), hold on,
                else
                    errorbar(1:size(Measurements,1),Measurements(:,counter_plot),zeros(size(Measurements,1),1),'gs','LineWidth',1,'MarkerSize',1, 'Color',[0.4 0.6 0]), hold on,
                end
            end
        end

        % Plot simulated data on top
        plot(1:size(Measurements,1),x(counter,:),'b.','MarkerSize',15)

        % Figure adjustment
        axis([0 size(Measurements,1)+1 0 1.1])
        set(gca,'XTick', [1:length(estim.Annotation)])
        set(gca,'XTickLabel', estim.Annotation)
        xtickangle(45)
        set(gca,'fontsize',15)
        t=title(state_names(counter));
        xt=xlabel('exp');
        y=ylabel('state-value');
        set(xt,'fontsize',15)
        set(y,'fontsize',15)
        set(t,'fontsize',15)
        hold off
        if ToSave
            saveas(h1,[Folder, '\', cell2mat(state_names(counter)) '_plot'],'tif')
        end

    end
end

if graphs(3) && sum(std(estim.Output_idx))==0
    % Plot optimal cost for each experiment
    hm=HeatMap(Diffs, 'RowLabels',1:size(estim.Output_idx,1),'ColumnLabels',estim.state_names(estim.Output_idx(1,:)),'Colormap',hot, 'Symmetric', false);
    addTitle(hm, 'Cross-error Analysis: Heatmap');
    if ToSave
        fighm=plot(hm);
        saveas(fighm,[Folder, '\CrossErrorHeatMap'],'tif');
        close(gcf)
    end
%%%% Modif removal of histogram
%     figure, hist(Diffs(:)); title('Cross-error Analysis: Histogram');
%     if ToSave
%         saveas(gcf,[Folder, '\CrossErrorHistogram'],'tif')
%     end
end

if graphs(4)
    h4=figure; hold on
    NrExps=length(estim.Output(:,1));
    for i=1:NrExps %for each exp
        subplot(NrExps,1,i)
        Y=state_names; X=Y;
        set(gca, 'xtick', 1:length(X))
        set(gca,'xticklabel',X)
        hold on
        Y2mean=x(:,i);
        Y2std=zeros(size(Y2mean));    
        errorbar(Y2mean, Y2std,'.b','LineWidth',1)
        hold on
        if ~isempty(SD)
            errorbar(Output_index(i,:), Measurements(i,:), SD(i,:),'.g', 'Color',[0.4 0.6 0])
        else
            errorbar(Output_index(i,:), Measurements(i,:), zeros(size(Measurements(i,:))),'.g', 'Color',[0.4 0.6 0])
        end
        axis([0, length(X)+1, 0, 1.1])
    end
    hold off
    if ToSave
        saveas(h4,[Folder, '\Measured_vs_Simul'],'tif')
    end
end

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


num_plots=length(Output_index);
NLines=ceil(sqrt(num_plots));
NCols=ceil(num_plots/NLines);

    h51=figure; hold on,
     suptitle('Dynamics through the simulation (outputs)');
    for p=1:size(Output_index,1)
        subplot(NLines, NCols,p)
        plot(squeeze(estim.AllofTheXs(p,1:T,Output_index(p,:))))
        axis([0, T+1, 0, 1]);
        hold off
    end
    legend(state_names(Output_index(p,:)))
       
    if ToSave
        saveas(h51,[Folder, '\ConvergenceOutputNodes'],'tif')
    end

    h52=figure; hold on,
      
        suptitle('Dynamics through the simulation (all nodes)')
    for p=1:size(Output_index,1)
        subplot(NLines, NCols,p)
        plot(squeeze(estim.AllofTheXs(p,1:T,:)))
        axis([0, T+1, 0, 1]);
        hold off
    end
      legend(state_names(:))
    if ToSave
        saveas(h52,[Folder, '\ConvergenceAllNodes'],'tif')
    end
end
MeanStateValueAll = x';
estim.Results.Optimisation.StdStateValueAll =  StdStateValueAll;
estim.Results.Optimisation.Diffs = Diffs;
end

