function FalconValidation(varargin)
% FalconValidation resimulates the optimized network and compare to the validation dataset
% FalconSimul(estim,bestx,ValidationDataset,FinalFolderName)

% :: Input values ::
% estim                     complete model definition
% bestx                     the vector of best optimised parameters
% ValidationDataset         name of the validation dataset
% FinalFolderName           name of the folder for saving results
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

estim=varargin{1}; optParams=varargin{2}; MeasFile=varargin{3};
ToSave=0;
if nargin>3, Folder=varargin{4}; ToSave=1; end

state_names=estim.state_names;

%=========================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the info from the measurement file

Point=find(ismember(MeasFile,'.'),1,'last'); %finding the last point in the file name
Ext=MeasFile(Point+1:end); %retrieving the extension
if strcmp(Ext,'txt') %if text file
    fid=fopen(MeasFile,'r');
    LineCounter=0;
    
    Input_vector=[];
    Input_index=[];
    Output_vector=[];
    Output_index=[];
    SD_vector=[];
    
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        
        LineCounter=LineCounter+1;
        ReadIO = regexp(tline,'\t','split');
        InputRaw=ReadIO(1);
        OutputRaw=ReadIO(2);
        ReadInput=strsplit(char(InputRaw),',');
        ReadOutput=strsplit(char(OutputRaw),',');
        
        if length(ReadIO)==3 % Shared index between output value and SD
            SDRaw=ReadIO(3);
            ReadSD=strsplit(char(SDRaw),',');
        end
        
        count_input=1;
        Input_idx_collect=[];
        Input_value_collect=[];
        for counter=1:(size(ReadInput,2)/2)
            idx_ReadInput=find(ismember(state_names,ReadInput(count_input)));
            value_ReadInput=str2num(cell2mat(ReadInput(count_input+1)));
            Input_idx_collect=[Input_idx_collect idx_ReadInput];
            Input_value_collect=[Input_value_collect value_ReadInput];
            count_input=count_input+2;
        end
        
        Input_vector=[Input_vector; Input_value_collect];
        Input_index=[Input_index; Input_idx_collect];
        
        count_output=1;
        Output_idx_collect=[];
        Output_value_collect=[];
        for counter=1:(size(ReadOutput,2)/2)
            idx_ReadOutput=find(ismember(state_names,ReadOutput(count_output)));
            value_ReadOutput=str2num(cell2mat(ReadOutput(count_output+1)));
            Output_idx_collect=[Output_idx_collect idx_ReadOutput];
            Output_value_collect=[Output_value_collect value_ReadOutput];
            count_output=count_output+2;
        end
        
        Output_vector=[Output_vector; Output_value_collect];
        Output_index=[Output_index; Output_idx_collect];
        
        if length(ReadIO)==3 % Shared index between output value and SD
            
            count_SD=1;
            SD_value_collect=[];
            for counter=1:(size(ReadSD,2)/2)
                value_ReadSD=str2num(cell2mat(ReadSD(count_SD+1)));
                SD_value_collect=[SD_value_collect value_ReadSD];
                count_SD=count_SD+2;
            end
            
            SD_vector=[SD_vector; SD_value_collect];
            
        end
        
    end
    fclose(fid);
elseif strcmp(Ext,'xls') || strcmp(Ext,'xlsx')
    [~,~,OtherIn]=xlsread(MeasFile,1);
    [~,~,OtherOut]=xlsread(MeasFile,2);
    [~,~,OtherErr]=xlsread(MeasFile,3);
    Input_index=[];
    Output_index=[];
    
     Input_vector=cell2mat(OtherIn(2:end,2:end)); 
    Annotation = (OtherIn(2:end,1)); 
    for jj=2:length(OtherIn(:,2)) %%%% modif Philippe OtherIn(:,1)
        Input_index_coll=[];
        for j=1:length(OtherIn(1,:))
            if ~isnan(cell2mat(OtherIn(jj,j)))
                Input_index_coll=[Input_index_coll,find(ismember(state_names,OtherIn(1,j)))];
            else
                Input_index_coll=[Input_index_coll,NaN];
            end
        end
        Input_index=[Input_index;Input_index_coll];
    end
    OtherOut(strcmp(OtherOut,'NaN'))={NaN}; %%%
    Output_vector=cell2mat(OtherOut(2:end,:));
    for jj=2:length(OtherOut(:,1))
        Output_index_coll=[];
        for j=1:length(OtherOut(1,:))
            if ~isnan(cell2mat(OtherOut(jj,j)))
                Output_index_coll=[Output_index_coll,find(ismember(state_names,OtherOut(1,j)))];
            else
                Output_index_coll=[Output_index_coll,NaN];
            end
        end
        Output_index=[Output_index;Output_index_coll];
    end
    OtherErr(strcmp(OtherErr,'NaN'))={NaN};
    if ~isempty(OtherErr)
        SD_vector=cell2mat(OtherErr(2:end,:));
    end
end

NewOutput_index=[];
Nexp=size(Output_index,1);
for oi=1:size(Output_index,2)
    thisId=Output_index(:,oi);
    thisId(isnan(thisId))=[];
    NewOutput_index=[NewOutput_index,unique(thisId)];
end
Output_index=repmat(NewOutput_index,Nexp,1);

%%% quick fix, search why SD_vector does not have the right length...
if ~isempty(SD_vector)
    SD_vector=SD_vector(:,1:size(Output_vector,2));
else
    SD_vector=zeros(size(Output_vector));
end
estim.SD=SD_vector;

% store the parameters in 'estim'
estim.Input=Input_vector;
estim.Input_idx=Input_index;
estim.Output=Output_vector;
estim.Output_idx=Output_index;

% =====================================

% Initialize necessary parameter and retrieve best parameter set (optParams)
tic
k=optParams;
StdStateValueAll=[];
MeanCostAll=[];
StdCostAll=[];
Diffs=[];
Max=10000;
% Input_index=estim.Input_idx;
% Output_index=estim.Output_idx;
Inputs=estim.Input;

Measurements=estim.Output;
SD=estim.SD;
param_index=estim.param_index;
ma=estim.ma;
mi=estim.mi;
N = numel(estim.Output)-sum(sum(isnan(estim.Output)));
np= numel(estim.param_vector);

% Perform simulation based on the best parameter set

n=estim.NrStates;
if isfield(estim,'MaxTime')
    MaxTime=estim.MaxTime;
else
    MaxTime=0;
end
%initial and successive number of steps for evaluation
if n<=25, initial_t=10; step_t=10;
elseif n>25 && n<=100, initial_t=100; step_t=10;
elseif n>100, initial_t=300; step_t=150;
end
TimeSoFar=tic;

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
x=rand(n,size(Measurements,1)); %initial random values for the nodes
x(Input_index(1,:),:)=Inputs';%(counter_exp,:); %fixing input nodes
xmeas=Measurements; %expected values for output nodes
break_point_ss=1;
runs=initial_t; % Initialize number of time steps

while break_point_ss
    if MaxTime>0
        if toc(TimeSoFar)>MaxTime
            break
        end
    end
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
mse=(sum(sum((xsim-xmeas).^2)))/N;
diff=mse;
disp(['MSE: ', num2str(mse), ' ; SSE: ', num2str(mse*N)])

x=x';
estim.Results.Optimisation.BestStates = x;
IdxOut=estim.Output_idx(1,:);
Heading=estim.state_names(IdxOut);
Numb=num2cell(x(:,IdxOut));
useexcel = isExcelPresent();
if useexcel
    xlswrite([Folder filesep 'Predictions.xls'],[Heading;Numb])
else
    
%     tab = table(x(:,IdxOut),'VariableNames',Heading);
%     writetable(tab,[FinalFolderName, filesep, 'Predictions.csv'],'Delimiter',',');
end

estim.Results.Validation.Outputs=Numb;


% PlotFitSummary: graph of state values at steady-state versus measurements in validation dataset
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
        errorbar(1:size(Measurements,1),Measurements(:,counter),zeros(size(Measurements,1),1),'gs','LineWidth',1,'MarkerSize',20/sqrt(num_plots),'Color',[0.4 0.6 0]), hold on,
    end

    % Plot simulated data on top
    plot(1:size(Measurements,1),x(:,Output_index(1,counter)),'b.','MarkerSize',20/sqrt(num_plots))

    % Figure adjustment
    axis([0 size(Measurements,1)+1 0 1.1])
    set(gca,'fontsize',15/sqrt(num_plots))
    set(gca,'XGrid','on')
    t=title(state_names(Output_index(1,counter)));
    xt=xlabel('experimental condition');
    set(xt,'fontsize',15/sqrt(num_plots))
    set(t,'fontsize',25/sqrt(num_plots))
    hold off
end

% scatter plot for observed vs estimated
h2=figure; hold on
for counter=1:num_plots
    subplot(NLines,NCols,counter), hold on,
    
    SSres=nansum((x(:,Output_index(1,counter))-Measurements(:,counter)).^2);
    SStot=nansum((Measurements(:,counter)-nanmean(Measurements(:,counter))).^2);

    plot(Measurements(:,counter), x(:,Output_index(1,counter)),'.k', 'MarkerSize', 10)
    plot([0,1],[0,1], 'Color', [0.8 0.8 0.8]);
    
    % Figure adjustment
    axis([0 1.1 0 1.1]), title([char(state_names(Output_index(1,counter))),': R^2= ', num2str(1-SSres/SStot)]);
    xlabel('observed');ylabel('estimated');
    grid on;
    hold off
end

h3=figure; hold on
X=estim.Output(:); Y=cell2mat(Numb(:));
Rsq2 = 1 - nansum((Y - X).^2)/nansum((Y - mean(Y)).^2);
plot(X,Y,'.k', 'MarkerSize', 10)
plot([0,1],[0,1], 'Color', [0.8 0.8 0.8]);

title(['Validation dataset: R^2=',num2str(Rsq2)])


if ToSave
    saveas(h1,[Folder,filesep,'Validation_plot'],'tif')
    saveas(h1,[Folder,filesep,'Validation_plot'],'fig')
    saveas(h2,[Folder,filesep,'Validation_plot2'],'tif')
    saveas(h2,[Folder,filesep,'Validation_plot2'],'fig')
    saveas(h3,[Folder,filesep,'Validation_plot3'],'tif')
    saveas(h3,[Folder,filesep,'Validation_plot3'],'fig')
    
end


end