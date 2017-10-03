function [ResultsLambda, ResultsParams] = Lcombined_analysis(estim,InputFile, MeasFileList, ContextsList, M, ListLambda12, ListLambda1group, PowerLambda12, PowerLambda1group, vartest ,FinalFolderName)

        
%% Analyis of combined L1 and L1/2 regularization
t= num2cell(vartest);
[PlotFitSummary, PlotFitIndividual, PlotHeatmapCost, PlotStateSummary, PlotStateEvolution] = deal(t{:});


%% Overview of the data analysis with general plots
% AIC evolution
figure, 
AICs=M(:,end-3);
AICs=reshape(AICs,length(ListLambda12),length(ListLambda1group));
imagesc(AICs);
title('AIC')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',num2str([0; PowerLambda1group(:)]))
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',num2str([0; PowerLambda12(:)]))
ylabel('L1 grouped')
xlabel('L1/2')

% BIC evolution
figure, 
BICs=M(:,end-2);
BICs=reshape(BICs,length(ListLambda12),length(ListLambda1group));
imagesc(BICs);
title('BIC')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',num2str([0, PowerLambda1group(:)]))
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',num2str([0, PowerLambda12(:)]))
ylabel('L1 grouped')
xlabel('L1/2')

% MSE evolution

figure,
MSEs=M(:,end-2);
MSEs=reshape(MSEs,length(ListLambda12),length(ListLambda1group));
imagesc(log(MSEs))
title('MSE L1 vs L1/2 regularization')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',['  0'; num2str(PowerLambda1group(:))])
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',['  0'; num2str(PowerLambda12(:))])
ylabel('L1 grouped')
xlabel('L1/2')

% parameter evolution
figure, 
Nps=M(:,end);
Nps=reshape(Nps,length(ListLambda12),length(ListLambda1group));
imagesc(Nps)
title('#Parameters')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',['  0'; num2str(PowerLambda1group(:))])
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',['  0'; num2str(PowerLambda12(:))])
ylabel('L1 grouped')
xlabel('L1/2')

%%% SBC
SBC= numel(estim.Output) .* log(M(:,1)) + M(:,end) .* log(numel(estim.Output));
figure, 
SBCs=reshape(SBC,length(ListLambda12),length(ListLambda1group));
imagesc(SBCs)
title('SBC')
set(gca,'XTick',1:length(ListLambda12))
set(gca,'xticklabels',num2str([0, PowerLambda1group(:)]))
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda1group))
set(gca,'yticklabels',num2str([0, PowerLambda12(:)]))
ylabel('L1 grouped')
xlabel('L1/2')


%% specific data analysis

% find minimal AIC
min_AIC = min(M(:,end-3));
AIC_pos = find(ismember(M(:,end-3),min_AIC));

% plot model fits for best AIC
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,M(AIC_pos,2:end-5),[PlotFitSummary, PlotFitIndividual, PlotHeatmapCost, PlotStateSummary, PlotStateEvolution],FinalFolderName);

best = M(AIC_pos,2:end-5);

best = best(estim.RegMatrix);
best_mean=mean(best, 2);
best_std=std(best, 0, 2);



[best_std, Order] = sortrows(best_std);
figure, 
bar(best_mean), hold on
errorbar(best_mean,best_std, '.k')
set(gca, 'xtick', 1:length(best_std))
set(gca, 'xticklabels', estim.param_vector(Order))
set(gca, 'xticklabelrotation', 90)
xlabel('parameters')
ylabel('parameter means')
title(['Parameter means and standard error at AIC=' num2str(min_AIC)])


%% identification of crucial parameters based on L1_groups
    
ref=FalconMakeModel(InputFile, MeasFileList{1}, 0.5)
Mgroup=M(1:length(ListLambda1group),2:end-5);

L=length(ContextsList);
N=size(Mgroup,2)/L;
STDs=zeros(size(Mgroup,1),N);
for g=1:size(Mgroup,1)
    for c=1:N
        Idx= 1+((c-1)*L):1+((c-1)*L)+L-1;
        STDs(g, c)=std(Mgroup(g, Idx));
    end
end

figure, imagesc(STDs)

Threshold=0.01;
STDs_pass=STDs<Threshold;
STDs_x=STDs_pass(1:end-1,:)~=STDs_pass(2:end,:);
[row col] = find(STDs_x);

ResultsLambda=[];
ResultsParams={};
for c=size(Mgroup,1):-1:1
    Idx=find(row==c);
    if ~isempty(Idx)
        ResultsLambda=[ResultsLambda; PowerLambda1group(c)];
        ResultsParams=[ResultsParams; {ref.param_vector(Idx)'}]
    end
end

end