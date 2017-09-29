function Lcombined_analysis(estim,ContextsList, M, ListLambda12, ListLambda1group, PowerLambda12, PowerLambda1group, vartest ,FinalFolderName)

%%%%%%%%%%%% TO DO %%%%%%%%%%%%
% % get the ranking of the most important parameters
% % additional selection criterion
%%%%%%%%%%%% %%%%%%%%%%%%


%% Analyis of combined L1 and L1/2 regularization
t= num2cell(vartest) 
[PlotFitSummary, PlotFitIndividual, PlotHeatmapCost, PlotStateSummary, PlotStateEvolution] = deal(t{:})


%% Overview of the data analysis with general plots
% AIC evolution
figure, 
AICs=M(:,end-3);
AICs=reshape(AICs,length(ListLambda12),length(ListLambda1group));
imagesc(AICs);
title('AIC')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',['0'; num2str(PowerLambda12(:))])
ylabel('L1 grouped')
xlabel('L1/2')

% BIC evolution
figure, 
BICs=M(:,end-2);
BICs=reshape(BICs,length(ListLambda12),length(ListLambda1group));
imagesc(BICs);
title('BIC')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',['0'; num2str(PowerLambda12(:))])
ylabel('L1 grouped')
xlabel('L1/2')

% MSE evolution

figure,
MSEs=M(:,end-2);
MSEs=reshape(MSEs,length(ListLambda12),length(ListLambda1group));
imagesc(log(MSEs))
title('MSE L1 vs L1/2 regularization')
ylabel('L1 grouped')
xlabel('L1/2')

% parameter evolution
figure, 
Nps=M(:,end);
Nps=reshape(Nps,length(ListLambda12),length(ListLambda1group));
imagesc(Nps)
title('#Parameters')
ylabel('L1 grouped')
xlabel('L1/2')

%%% SBC
SBC= numel(estim.Output) .* log(M(:,1)) + M(:,end) .* log(numel(estim.Output));
figure, 
SBCs=reshape(SBC,length(ListLambda12),length(ListLambda1group));
imagesc(SBCs)
title('SBC')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',['0'; num2str(PowerLambda12(:))])
ylabel('L1 grouped')
xlabel('L1/2')


%% specific data analysis

% find minimal AIC
min_AIC = min(M(:,end-3));
AIC_pos= find(ismember(M(:,end-3),min_AIC));

% plot model fits for best AIC
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,M(AIC_pos,2:end-5),[PlotFitSummary, PlotFitIndividual, PlotHeatmapCost, PlotStateSummary, PlotStateEvolution],FinalFolderName);


 
best_dataLambda = M(AIC_pos,2:end-5);

k=1;
counter=1;
parameterValues_mean=[];
parameterValues_std=[];
for k=1:length(best_dataLambda)/length(ContextsList)
    data=best_dataLambda(1,counter:counter+length(ContextsList)-1);
    test_mean_std(k,1)=mean(data);
    test_mean_std(k,2)=std(data);
    counter=counter+length(ContextsList);
    k=k+1;
end
test_mean_std=sortrows(test_mean_std);
figure; errorbar(test_mean_std(:,1),test_mean_std(:,2))
xlabel('parameters')
ylabel('parameter means')
title(['Parameter means and standard error at AIC=' num2str(min_AIC)])  


%% identification of crucial parameters based on L1_groups

% find minimal AIC
min_AIC = min(M(:,end-3));
AIC_pos= find(ismember(M(:,end-3),min_AIC));
 

Full_M = M(:,2:end-5);
%dataLambda4u16=M(410,2:196)

parameterValues_mean=[];
parameterValues_std=[];
for rows= 1: size(M,1)
    k=1;
    counter=1;
    for k=1:length(Full_M(rows,:))/length(ContextsList)
        data=Full_M(rows,counter:counter+length(ContextsList)-1);
        test_mean(rows,k)=mean(data);
        test_std(rows,k)=std(data);
        counter=counter+length(ContextsList);
        k=k+1;
    end
end
% test_std_sort=sortrows(test_std(:,1))

figure; plot(test_std(1:length(ListLambda1group),:));
ylabel('standard deviation')
xlabel('Lambda')
title(['Regularization dynamics'])  

% get for each param the value where STD  is below 1%

[row col] = find(test_std(1:length(ListLambda1group),:)< 0.01)
A= [row col]


end