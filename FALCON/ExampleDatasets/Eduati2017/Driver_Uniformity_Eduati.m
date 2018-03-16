
%%
%%% REG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, clc, %close all force

% Define optmisation options
optRound=4; % Number of optimisation round
MaxFunEvals=50000; % Number of maximal function being evaluated (3000 = default)
MaxIter=50000; % Number of maximal iteration being evaluated (3000 = default)
Parallelisation=1; % Use multiple cores for optimisation? (0=no, 1=yes)
HLbound=0.5; % Qualitative threshold between high and low inputs
InitIC=2; % Initialise parameters' distribution (1=uniform, 2=normal)

% Define plotting and saving (0=no, 1=yes)
PlotFitEvolution    = 0; % Graph of optimise fitting cost over iteration
PlotFitSummary      = 1; % Graph of state values at steady-state versus measurements (all in 1)
PlotFitIndividual   = 0; % Graph of state values at steady-state versus measurements (individual)
PlotHeatmapCost     = 1; % Heatmaps of optimal costs for each output for each condition absolute cost
PlotStateSummary    = 0; % Graph of only state values at steady-sate (all in 1)
PlotStateEvolution  = 0; % Graph of state values evolution over the course of the simulation (two graphs)
PlotBiograph        = 0; % Graph of network topology, nodes activities, and optimised parameters
PlotAllBiographs    = 0; % (Only for machines with strong GPUs) Plot all Biographs above

N=0;
ExpName='test_'

FALCONFolder = fileparts(which('DriverFalcon'));
ExampleFolder = [FALCONFolder filesep 'ExampleDatasets' filesep 'Eduati2017'];
InputFile=[ExampleFolder filesep 'Eduati_net2.xlsx']; %% Network structure file
ContextsList={'CAR1','CCK81','COLO320HSR', 'DIFI', 'HCT116', 'HT29', ...
    'HT115', 'SKCO1', 'SNUC2B', 'SNUC5', 'SW620', 'SW837',...
    'SW1116', 'SW1463'};
MeasFileList={};
for f=1:length(ContextsList)
      MeasFileList=[MeasFileList,[ExampleFolder filesep char(ContextsList(f)) '_noCTRL_meas.xlsx']]; 
end
GlobalEdgesList=[ExampleFolder filesep 'Eduati_nofixed.xlsx']; 
% Build a Global FALCON model for optimisation

[estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);

NumFixedInteractions=0;

%%%%
%%%%

% Build a Global FALCON model for optimisation
[estim, File, GlobalFileName] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
% Define optimisation options
estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
estim.SSthresh=0.001;
Np=length(estim.param_vector(NumFixedInteractions+1:end)); Nc=numel(ContextsList); Ns=length(estim.state_names);
R=[];

PowerLambdaCluster=-20:1:4;

%%%%
%%%%

ListLambdaCluster=2.^PowerLambdaCluster;

AICs=zeros(length(ListLambdaCluster),1);
MSEs=zeros(length(ListLambdaCluster),1);
BICs=zeros(length(ListLambdaCluster),1);
NPs=zeros(length(ListLambdaCluster),1);
Times=zeros(length(ListLambdaCluster),1);
%%%%%%%%%%%%%%
Params=zeros(length(ListLambdaCluster),Np);
StateValues=zeros(length(ListLambdaCluster),Ns);

for c=1:(Nc):Np
    R=[R;reshape(1:(Nc),1,Nc)+(c-1)];

end

estim.Reg='LCluster';
estim.RegMatrix.Cluster=R+NumFixedInteractions;

%%%%%%

l1=0;
Dims=[];
for L1=ListLambdaCluster
    l1=l1+1;
    N=N+1;
    if AICs(l1)==0
        Dims=[Dims;L1];
        % Create a save folder
        SaveFolderName=['Results', filesep, ExpName, '_', char(datetime)];
        FinalFolderName=strrep(SaveFolderName, ':', '.');
        mkdir(FinalFolderName) % Automatically generate a folder for saving

        estim.Lambda=L1;

        fxt_all=[];
        Costs_all=[];
        parfor cc=1:optRound
            tic
            k=FalconIC(estim,'normal'); %initial conditions
            [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
            Fitting_Time=toc;

            fxt_all=[fxt_all;[fval xval Fitting_Time]];
            Costs_all=[Costs_all;[AIC, MSE, BIC,Np]];
            beep; pause(0.5); beep;
        end

        [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
        [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

        Idx=find(ismember(min(Costs_all(:,end-2)),Costs_all(:,end-2)));

        AICs(l1)=Costs_all(Idx,end-3);
        MSEs(l1)=Costs_all(Idx,end-2);
        BICs(l1)=Costs_all(Idx,end-1);
        NPs(l1)=Costs_all(Idx,end);
        Times(l1)=fxt_all(Idx,end);
        Params(l1,:)=bestx;
        States(l1,:)=MeanStateValueAll(:);

        save(ExpName)
    end
end

bestx(estim.RegMatrix.Cluster)

figure,
subplot(2,2,1)
plot(AICs,'k', 'linewidth', 2), set(gca, 'xtick', 1:2:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster(1:2:end)}), set(gca, 'xticklabelrotation', 90), title('AIC')
xlabel('log_2(\lambda_{cluster})')
subplot(2,2,2)
plot(MSEs,'k', 'linewidth', 2), set(gca, 'xtick', 1:2:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster(1:2:end)}), set(gca, 'xticklabelrotation', 90), title('MSE')
xlabel('log_2(\lambda_{cluster})')
subplot(2,2,3)
plot(BICs,'k', 'linewidth', 2), set(gca, 'xtick', 1:2:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster(1:2:end)}), set(gca, 'xticklabelrotation', 90), title('BIC')
xlabel('log_2(\lambda_{cluster})')
subplot(2,2,4)
plot(NPs,'k', 'linewidth', 2), set(gca, 'xtick', 1:2:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster(1:2:end)}), set(gca, 'xticklabelrotation', 90), title('#Parameters')
xlabel('log_2(\lambda_{cluster})')
suptitle('Model characteristics across regularization strengths')

M=estim.RegMatrix.Cluster;
figure,
Colors=distinguishable_colors(size(M,1));
for group=1:size(M,1)
    P=Params(:,M(group,:));
    plot(1:size(P,1),P,'Color', Colors(group,:), 'Linewidth', 2), hold on
end
set(gca,'xtick',1:size(Params,1))
set(gca,'xticklabel', {PowerLambdaCluster})
set(gca, 'xticklabelrotation', 90)
title('Regularization path')

save(ExpName)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clear all, load('20180117_reg')

IdxBest=find(AICs==min(AICs))

% Create a save folder
SaveFolderName=['FinalResults_', ExpName, '_', char(datetime)];
FinalFolderName=strrep(SaveFolderName, ':', '.');
mkdir(FinalFolderName) % Automatically generate a folder for saving

bestx=Params(IdxBest,:)
bestx(M)

estim=FalconMakeModel('Expanded_Corrected.txt',GlobalFileName,HLbound)
estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
estim.SSthresh=0.001;

fxt_all=[];
Costs_all=[];
parfor cc=1:optRound
    tic
    k=FalconIC(estim,'normal'); %initial conditions
    [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
    Fitting_Time=toc;

    fxt_all=[fxt_all;[fval xval Fitting_Time]];
    Costs_all=[Costs_all;[AIC, MSE, BIC,Np]];
    beep; pause(0.5); beep;
end

[bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

InferredParams=([bestx([1 3 2 4 5 6 7 8 10]); bestx([ 1 3 2 4 5 6 7 8 10]); bestx([ 1 3 2 4 5 6 7 9 11]); bestx([ 1 3 2 4 5 6 7 9 11])])

Measured=estim.Output; Simulated=MeanStateValueAll(:,estim.Output_idx(1,:));

CostsPerLine=[];
for c=0:3
    n=2*c;
    Num=numel(Measured(:,1+n:2+n))
    thisCost=(Measured(:,1+n:2+n) - Simulated(:,1+n:2+n)).^2;
    CostsPerLine(c+1)=mean(thisCost(:));
    
end
figure, 
subplot(2,2,1),
bar(CostsPerLine, 'FaceColor', [0.8 0.8 0.8])
title('Fitting cost (MSE)')
set(gca, 'xtick', 1:4)
set(gca, 'xticklabel', {'CL1','CL2','CL3','CL4'})
subplot(2,2,2),
plot(Measured, Simulated,'.k', 'MarkerFaceColor', [0.8 0.8 0.8], 'Markersize', 8), hold on
xlabel('Measured')
ylabel('Simulated')
X=Measured(:);
X=[ones(length(X),1) X];
b=X\Simulated(:);
LM=fitlm(Measured(:), Simulated(:));
rho=LM.Rsquared.Adjusted;
plot([0; 1], b, '-.k')
axis([0 1.1 0 1.1])
title(['Data correlation plot, r^2_{adj}=', num2str(rho)])
subplot(2,2,3),
plot(RealParams(:), InferredParams(:),'.k', 'MarkerFaceColor', [0.8 0.8 0.8], 'Markersize', 8), hold on
xlabel('Actual')
ylabel('Estimated')
X=InferredParams(:);
X=[ones(length(X),1) X];
b=X\RealParams(:);
LM=fitlm(RealParams(:), InferredParams(:));
rho=LM.Rsquared.Adjusted;
plot([0; 1], b, '-.k')
axis([0 1.1 0 1.1])
title(['Parameter correlation plot, r^2_{adj}=', num2str(rho)])

save('Toy_0117')


%%





















%%
%%% 2D-REG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, clc,

% Define optmisation options
optRound=4; % Number of optimisation round
MaxFunEvals=50000; % Number of maximal function being evaluated (3000 = default)
MaxIter=50000; % Number of maximal iteration being evaluated (3000 = default)
Parallelisation=1; % Use multiple cores for optimisation? (0=no, 1=yes)
HLbound=0.5; % Qualitative threshold between high and low inputs
InitIC=2; % Initialise parameters' distribution (1=uniform, 2=normal)

% Define plotting and saving (0=no, 1=yes)
PlotFitEvolution    = 0; % Graph of optimise fitting cost over iteration
PlotFitSummary      = 0; % Graph of state values at steady-state versus measurements (all in 1)
PlotFitIndividual   = 0; % Graph of state values at steady-state versus measurements (individual)
PlotHeatmapCost     = 0; % Heatmaps of optimal costs for each output for each condition absolute cost
PlotStateSummary    = 0; % Graph of only state values at steady-sate (all in 1)
PlotStateEvolution  = 0; % Graph of state values evolution over the course of the simulation (two graphs)
PlotBiograph        = 0; % Graph of network topology, nodes activities, and optimised parameters
PlotAllBiographs    = 0; % (Only for machines with strong GPUs) Plot all Biographs above

N=0;

% fix Lambda values
PowerLambda12=-20:1:-1; % L1/2
% PowerLambda1group=-6:0.5:-1; % L1group
% PowerLambda1smooth=-6:0.5:-1; % L1smooth

PowerLambdaCluster=-20:1:-1;


% List12=10.^PowerLambda12;
% Listgroup=10.^PowerLambda1group;
% Listsmooth=10.^PowerLambda1smooth;
ListLambda12=2.^PowerLambda12;
ListLambdaCluster=2.^PowerLambdaCluster;

AICs=zeros(length(ListLambda12),length(ListLambdaCluster));
MSEs=zeros(length(ListLambda12),length(ListLambdaCluster));
BICs=zeros(length(ListLambda12),length(ListLambdaCluster));
NPs=zeros(length(ListLambda12),length(ListLambdaCluster));
Times=zeros(length(ListLambda12),length(ListLambdaCluster));
Params=zeros(length(ListLambda12),length(ListLambdaCluster),32);

CellLines=1:4;
NumFixedInteractions=0;

MeasFileList={};
ContextsList={};
for c=CellLines

    MeasFileList=[MeasFileList,['ToyCluster_Data' num2str(c) '.xlsx']];
    ContextsList=[ContextsList, ['C' num2str(c)]];

end

InputFile='Toy_net3.xlsx'; %% Network structure file
ExpName='Toy1214_tryReg_2d';
GlobalEdgesList='Toy_net1_empty.xlsx';

% Build a Global FALCON model for optimisation
[estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
% Define optimisation options
estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
estim.SSthresh=0.001;
Np=length(estim.param_vector(NumFixedInteractions+1:end)); Nc=numel(CellLines);
R=[];
for c=1:(Nc):Np
    R=[R;reshape(1:(Nc),1,Nc)+(c-1)];

end

estim.Reg='PruneCluster';
estim.RegMatrix.Cluster=R+NumFixedInteractions;

%%%%%%

l1=0;
Dims=[];
for L1=ListLambda12
    l2=0;
    l1=l1+1;
    for L2=ListLambdaCluster
        l2=l2+1;
        N=N+1;
        if AICs(l1,l2)==0
            Dims=[Dims;[L1,L2]];
            % Create a save folder
            SaveFolderName=['Results_', ExpName, '_', char(datetime)];
            FinalFolderName=strrep(SaveFolderName, ':', '.');
            mkdir(FinalFolderName) % Automatically generate a folder for saving

            estim.Lambda=[L1,L2];

            fxt_all=[];
            Costs_all=[];
            parfor cc=1:optRound
                tic
                k=FalconIC(estim,'normal'); %initial conditions
                [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
                Fitting_Time=toc;

                fxt_all=[fxt_all;[fval xval Fitting_Time]];
                Costs_all=[Costs_all;[AIC, MSE, BIC,Np]];
                beep; pause(0.5); beep;
            end

            [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
            [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

            Idx=find(ismember(min(Costs_all(:,end-2)),Costs_all(:,end-2)));

            AICs(l1,l2)=Costs_all(Idx,end-3);
            MSEs(l1,l2)=Costs_all(Idx,end-2);
            BICs(l1,l2)=Costs_all(Idx,end-1);
            NPs(l1,l2)=Costs_all(Idx,end);
            Times(l1,l2)=fxt_all(Idx,end);
            Params(l1,l2,:)=bestx;

            save('Toy_1214_cluster2d_addedge')
        end
    end
end

bestx(estim.RegMatrix.Cluster)

figure,
subplot(2,2,1)
imagesc(AICs), set(gca, 'xtick', 1:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster}), set(gca, 'xticklabelrotation', 90), title('AIC')
set(gca, 'ytick', 1:length(ListLambda12)), set(gca, 'yticklabel', {PowerLambda12}), colorbar, colormap('hot')
subplot(2,2,2)
imagesc(log(MSEs)), set(gca, 'xtick', 1:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster}), set(gca, 'xticklabelrotation', 90), title('log(MSE)')
set(gca, 'ytick', 1:length(ListLambda12)), set(gca, 'yticklabel', {PowerLambda12}), colorbar, colormap('hot')
subplot(2,2,3)
imagesc(BICs), set(gca, 'xtick', 1:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster}), set(gca, 'xticklabelrotation', 90), title('BIC')
set(gca, 'ytick', 1:length(ListLambda12)), set(gca, 'yticklabel', {PowerLambda12}), colorbar, colormap('hot')
subplot(2,2,4)
imagesc(NPs), set(gca, 'xtick', 1:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster}), set(gca, 'xticklabelrotation', 90), title('Number of Params')
set(gca, 'ytick', 1:length(ListLambda12)), set(gca, 'yticklabel', {PowerLambda12}), colorbar, colormap('hot')
figure, 
imagesc(Times), set(gca, 'xtick', 1:length(ListLambdaCluster)), set(gca, 'xticklabel', {PowerLambdaCluster}), set(gca, 'xticklabelrotation', 90), title('CPU time')
set(gca, 'ytick', 1:length(ListLambda12)), set(gca, 'yticklabel', {PowerLambda12}), colorbar, colormap('hot')





M=estim.RegMatrix.Cluster;
mn=length(ListLambda12);
m=floor(sqrt(mn)); n=ceil(mn/m);
figure,
for c=1:length(ListLambda12)
    subplot(m,n,c)
    Colors=distinguishable_colors(size(M,1));
    for group=1:size(M,1)
        P=squeeze(Params(c,:,M(group,:)));
        plot(1:size(P,1),P,'Color', Colors(group,:), 'Linewidth', 1), hold on
    end
    set(gca,'xtick',1:size(Params,1))
    set(gca,'xticklabel', {PowerLambdaCluster})
    set(gca, 'xticklabelrotation', 90)
    title(['\lambda_{1/2}= 2^{', num2str(PowerLambda12(c)), '}'])
end
suptitle('Regularization paths')




ParamDiffs=Params(:,2:end,:)-Params(:,1:end-1,:);
figure,
for c=1:length(ListLambda12)
    subplot(m,n,c)
    Colors=distinguishable_colors(size(M,1));
    for group=1:size(M,1)
        P=squeeze(ParamDiffs(c,:,M(group,:)));
        plot(1:size(P,1),P,'Color', Colors(group,:), 'Linewidth', 1), hold on
    end
    ylim([-1 1])
    set(gca,'xtick',1:size(ParamDiffs,1))
    set(gca,'xticklabel', {PowerLambdaCluster(2:end)})
    set(gca, 'xticklabelrotation', 90)
    title(['\lambda_{1/2}= 2^{', num2str(PowerLambda12(c)), '}'])
end
suptitle('Regularization paths, diff versus unregularized')




%%
















%%% visualize
%%% AICs
figure, imagesc(squeeze(AICs(:,:,1))), ylabel('L12'), xlabel('Lgroup'), 
colorbar, colormap('hot'), title('AIC, time-variant model')
set(gca, 'xtick', 1:length(Listgroup)), set(gca, 'xticklabel', Listgroup), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(AICs(:,1,:))), ylabel('L12'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('AIC, context-independent')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(AICs(1,:,:))), ylabel('Lgroup'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('AIC, unpruned topology')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(Listgroup)), set(gca, 'yticklabel', Listgroup)

Z=permute(AICs, [3 2 1]); Z=Z(:); Sz=-Z; Sz=(Sz-min(Sz))/(max(Sz)-min(Sz))+eps;
Dims(Dims==0)=10^-7
LogDims=log10(Dims)
figure, scatter3(LogDims(:,1), LogDims(:,2), LogDims(:,3),Sz*100,Z,'filled')
xlabel('log(Lambda_{prune})')
ylabel('log(Lambda_{group})')
zlabel('log(Lambda_{smooth})')
colormap('hot')
cb = colorbar;
cb.Label.String = 'AIC';

figure, hist(AICs(:),50), title('Histogram of AIC values')

Best_AIC=min(AICs(:));
[xslice,yslice,zslice]=ind2sub(size(AICs),find(AICs == Best_AIC));

figure,
slice(1:length(List12),1:length(Listgroup),1:length(Listsmooth),AICs,xslice,yslice,zslice)
xlabel('Pruning')
ylabel('Grouping')
zlabel('Smoothing')

colormap hot









%%% MSEs
figure, imagesc(squeeze(MSEs(:,:,1))), ylabel('L12'), xlabel('Lgroup')
colorbar, colormap('hot'), title('MSE, time-variant model')
set(gca, 'xtick', 1:length(Listgroup)), set(gca, 'xticklabel', Listgroup), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(MSEs(:,1,:))), ylabel('L12'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('MSE, context-independent')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(MSEs(1,:,:))), ylabel('Lgroup'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('MSE, unpruned topology')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(Listgroup)), set(gca, 'yticklabel', Listgroup)

%%% BICs
figure, imagesc(squeeze(BICs(:,:,1))), ylabel('L12'), xlabel('Lgroup')
colorbar, colormap('hot'), title('BIC, time-variant model')
set(gca, 'xtick', 1:length(Listgroup)), set(gca, 'xticklabel', Listgroup), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(BICs(:,1,:))), ylabel('L12'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('BIC, context-independent')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(BICs(1,:,:))), ylabel('Lgroup'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('BIC, unpruned topology')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(Listgroup)), set(gca, 'yticklabel', Listgroup)

%%% NPs
figure, imagesc(squeeze(NPs(:,:,1))), ylabel('L12'), xlabel('Lgroup')
colorbar, colormap('hot'), title('Number of parameters, time-variant model')
set(gca, 'xtick', 1:length(Listgroup)), set(gca, 'xticklabel', Listgroup), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(NPs(:,1,:))), ylabel('L12'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('Number of parameters, context-independent')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(NPs(1,:,:))), ylabel('Lgroup'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('Number of parameters, unpruned topology')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(Listgroup)), set(gca, 'yticklabel', Listgroup)

%%% Times
figure, imagesc(squeeze(Times(:,:,1))), ylabel('L12'), xlabel('Lgroup')
colorbar, colormap('hot'), title('CPU time, time-variant model')
set(gca, 'xtick', 1:length(Listgroup)), set(gca, 'xticklabel', Listgroup), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(Times(:,1,:))), ylabel('L12'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('CPU time, context-independent')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(List12)), set(gca, 'yticklabel', List12)
figure, imagesc(squeeze(Times(1,:,:))), ylabel('Lgroup'), xlabel('Lsmooth')
colorbar, colormap('hot'), title('CPU time, unpruned topology')
set(gca, 'xtick', 1:length(Listsmooth)), set(gca, 'xticklabel', Listsmooth), 
set(gca, 'ytick', 1:length(Listgroup)), set(gca, 'yticklabel', Listgroup)








%%%%%%%%%%%
%%
load('ResultsFinallyHigh.mat')
B3=bestx
load('ResultsFinally.mat')
B2=bestx
load('ResultsFinallyControl.mat')
B=bestx
BB=[B', B2', B3']
BB=BB(6:end,:)

%%
%%%across timepoints
for t=1:3:size(BB,1)
    figure, hold on
    data=BB(t:t+2,:)
    plot(data')
    title((estim.param_vector(t+5)))
    legend({'1h','24h','48h'})
    xlabel('regularization strength')
    set(gca, 'xtick', 1:3)
    set(gca, 'xticklabel', {'None','Low','High'} )
    axis([0.9 3.1 0 1.05])
end

%%%across cell lines
for t=1:18:size(BB,1)
    for t2=0:2
        figure, hold on
        data=BB((t+t2):3:(t+t2)+18,:)
        plot(data')
        title((estim.param_vector(t)))
        legend({'cl1','cl2','cl3', 'cl4', 'cl5', 'cl6'})
        xlabel('regularization strength')
        set(gca, 'xtick', 1:3)
        set(gca, 'xticklabel', {'None','Low','High'} )
        axis([0.9 3.1 0 1.05])
    end
    
    figure, hold on
    data=BB(t:t+17,:)
    plot(data')
    title((estim.param_vector(t)))
    xlabel('regularization strength')
    set(gca, 'xtick', 1:3)
    set(gca, 'xticklabel', {'None','Low','High'} )
    axis([0.9 3.1 0 1.05])
end





%%
% MeasFileList={'Sorger_Data_mod_All.xlsx'};
% ContextsList='ALL';

Max=10;
GlobalEdgesList='SorgerNet_fixed.xlsx';


%%% no regularization


    
% Create a save folder
SaveFolderName=['Results_', ExpName, '_noReg_', char(datetime)];
FinalFolderName=strrep(SaveFolderName, ':', '.');
mkdir(FinalFolderName) % Automatically generate a folder for saving

% Build a Global FALCON model for optimisation
[estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);

% Define optimisation options
estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
estim.SSthresh=0.001;

estim.Lambda=2^-6;
estim.Reg='L1/2';
estim.Reg='none'
fxt_all=[];
Costs_all=[];
parfor cc=1:optRound
    tic
    k=FalconIC(estim,'normal'); %initial conditions
    [xval,fval,AIC,MSE, BIC,Np]=FalconObjFun(estim,k); %objective function
    Fitting_Time=toc;

    fxt_all=[fxt_all;[fval xval Fitting_Time]];
    Costs_all=[Costs_all;[AIC, MSE, BIC,Np]];
    beep; pause(0.5); beep;
end

[bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
[MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);

X=estim.Output(:); Y=MeanStateValueAll(:,estim.Output_idx(1,:)); Y=Y(:);

s=size(estim.Output,1);
figure, plot(X,Y,'.k'), hold on, plot(X(s*2+1:s*4),Y(s*2+1:s*4),'.r'), hold on, plot(X(s*4+1:s*6),Y(s*4+1:s*6),'.g')
save('L-Noreg')


%% L1/2 regularization

if L1_2
    
    M=[];
    
    for L=List12
        %     Create a save folder
        SaveFolderName=['Results_', ExpName, '_L12reg_', char(datetime)];
        FinalFolderName=strrep(SaveFolderName, ':', '.');
        mkdir(FinalFolderName) % Automatically generate a folder for saving
        
%         Build a Global FALCON model for optimisation
        [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
        
%         Define optimisation options
        estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
        estim.SSthresh=0.001;
        
        estim.Lambda=L;
        estim.Reg='L1/2';
        fxt_all=[];
        Costs_all=[];
        for cc=1:optRound
            tic
            k=FalconIC(estim,'normal'); %initial conditions
            [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
            Fitting_Time=toc;
            
            fxt_all=[fxt_all;[fval xval Fitting_Time]];
            Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
            beep; pause(0.5); beep;
        end
        
        [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
        [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
        idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
        M=[M;fxt_all(idx,:),Costs_all(idx,:)];
    end
    
    figure,
    subplot(3,1,1)
    plot(M(:,2:end-5)- repmat(M(1,2:end-5),size(M,1),1))
    set(gca,'XTick',1:length(List12))
    set(gca,'xticklabels',['0'; num2str(PowerLambda12(:))])
    set(gca,'xticklabelrotation',90)
    subplot(3,1,2)
    plot(M(:,[2:5,34:37]))
    set(gca,'XTick',1:length(List12))
    set(gca,'xticklabels',['0'; num2str(PowerLambda12(:))])
    set(gca,'xticklabelrotation',90)
    subplot(3,1,3)
    bar(M(:,end-5))
    set(gca,'XTick',1:length(List12))
    set(gca,'xticklabels',['0'; num2str(PowerLambda12(:))])
    set(gca,'xticklabelrotation',90)
    legend('AIC')
    
    
    
end

save('L12Large')


%% L1group regularization
if L1Groups
    
    M=[];
    
    for L=Listgroup
        %     Create a save folder
        SaveFolderName=['Results_', ExpName, '_L1Groupreg_', char(datetime)];
        FinalFolderName=strrep(SaveFolderName, ':', '.');
        mkdir(FinalFolderName) % Automatically generate a folder for saving
        
        %     Build a Global FALCON model for optimisation
        [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
        
        %     Define optimisation options
        estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
        estim.SSthresh=0.001;
        
        estim.Lambda=L;
        estim.RegMatrix=(reshape(1:length(estim.param_vector),length(ContextsList), length(estim.param_vector)/length(ContextsList)))';
        estim.Reg='L1Groups';
        fxt_all=[];
        Costs_all=[];
        for cc=1:optRound
            tic
            k=FalconIC(estim,'normal'); %initial conditions
            [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
            Fitting_Time=toc;
            
            fxt_all=[fxt_all;[fval xval Fitting_Time]];
            Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
            beep; pause(0.5); beep;
        end
        
        [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
        [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
        idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
        M=[M;fxt_all(idx,:),Costs_all(idx,:)];
    end
    
    figure,
    subplot(4,1,1)
    plot(M(:,2:end-5)-repmat(M(1,2:end-5), size(M,1),1))
    set(gca,'XTick',1:length(Listgroup))
    set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
    set(gca,'xticklabelrotation',90)
    title(' All Parameters')
    
    subplot(4,1,2)
    P=(M(:,2:end-5))';
    P2=(std(reshape(P(:),4,length(P(:))/length(ContextsList))))';
    S=reshape(P2, length(P2)/size(P,2), size(P,2));
    Collapsed=sum(S<0.01);
    Nparams=(max(Collapsed)-Collapsed).*length(ContextsList) + Collapsed;
    bar(Nparams)
    set(gca,'XTick',1:length(Listgroup))
    set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
    set(gca,'xticklabelrotation',90)
    title('#Parameters')
    
    subplot(4,1,3)
    plot(M(:,6:9));
    set(gca,'XTick',1:length(Listgroup))
    set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
    set(gca,'xticklabelrotation',90)
    title('Example parameter k2')
    legend({'Cell Line a','Cell Line b','Cell Line c','Cell Line d'})
    
    subplot(4,1,4)
    bar(M(:,end-5))
    set(gca,'XTick',1:length(Listgroup))
    set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
    set(gca,'xticklabelrotation',90)
    title('AIC')
    
    figure, imagesc(P)
    set(gca,'XTick',1:length(Listgroup))
    set(gca,'xticklabels',['0'; num2str(PowerLambda1group(:))])
    set(gca,'xticklabelrotation',90)
      
    save('LGroups')
    
end



%% Ldrug regularizaiton

if Ldrug
    
    
    M=[];
    count= 0
    for L1=List12
        count = count + 1
        for L2=Listgroup
            % Create a save folder
            SaveFolderName=['Results_', ExpName, '_L1Drugreg_', char(datetime)];
            FinalFolderName=strrep(SaveFolderName, ':', '.');
            mkdir(FinalFolderName) % Automatically generate a folder for saving
            
            % Build a Global FALCON model for optimisation
            [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
            % Define optimisation options
            estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
            estim.SSthresh=0.001;
            
            estim.Lambda=[L1,L2];
            estim.RegMatrix=(reshape(1:length(estim.param_vector),length(ContextsList), length(estim.param_vector)/length(ContextsList)))';
            estim.Reg='Ldrug';
            fxt_all=[];
            Costs_all=[];
            parfor cc=1:optRound
                tic
                k=FalconIC(estim,'normal'); %initial conditions
                [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
                fxt_all=[fxt_all;[fval xval toc]];
                Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
                beep; pause(0.5); beep;
            end
            
            [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
            [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
            idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
            M=[M;fxt_all(idx,:),Costs_all(idx,:)];
        end
    end
    
    save('LCombine')
    
    Lcombined_analysis(estim, ContextsList, M, List12, Listgroup, PowerLambda12, PowerLambda1group, [PlotFitSummary, PlotFitIndividual, PlotHeatmapCost, PlotStateSummary, PlotStateEvolution], FinalFolderName)
    
end
