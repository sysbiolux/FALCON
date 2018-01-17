% ======================================
% Driver script (run-through) for FALCON
% ======================================

% FalconInstall % In case the Falcon toolbox has not yet been added to Matlab's path
clc, clear all % clear screen and workspace

% Define optmisation options
optRound=4; % Number of optimisation round
MaxFunEvals=8000; % Number of maximal function being evaluated (3000 = default)
MaxIter=8000; % Number of maximal iteration being evaluated (3000 = default)
% Parallelisation=1; % Use multiple cores for optimisation? (0=no, 1=yes)
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
% Additional analyses after the optimisation with the default setting (0=no, 1=yes)
Resampling_Analysis = 0; % Resampling of experimental data and re-optimise
NDatasets           = 10;% Number of artificial datasets from which to resample.
LPSA_Analysis       = 0; % Local parameter sensitivity analysis
Fast_Option         = 1; % Performing faster LPSA by stopping if fitting costs go over a set threshold value
LPSA_Increments     = 10; % Number of increments for LPSA. Increase for finer resolution
KO_Analysis         = 0; % Parameter knock-out analysis

% Read model and measurement files

PowerLambda=-10:0.5:10;
ListLambda=[0,10.^PowerLambda];
% ExampleFolder = ['Eduati' filesep 'DataFiles'];
% InputFile= [ExampleFolder filesep 'Eduati_net2.xlsx'];

FALCONFolder = fileparts(which('DriverFalcon'));
ExampleFolder = [FALCONFolder filesep 'ExampleDatasets' filesep 'Greta'];
InputFile=[ExampleFolder filesep 'Nicole_gene_inhibitor_network_gating_FOXO_big.xlsx']; %% Network structure file
ExpName='Nicole_InhibitorTreatment_Network';
ContextsList={'ctrl','MTX','PDT'};MeasFileList={};
for f=1:length(ContextsList)
      MeasFileList=[MeasFileList,[ExampleFolder filesep 'Nicole_crosstalk_' char(ContextsList(f)) '_inhib_big.xlsx']]; 
end
Max=10; 
GlobalEdgesList=[ExampleFolder filesep 'Nicole_fixed.xlsx']; 


%%
% 
% % Create a save folder
% SaveFolderName=['Results_', ExpName, '_noReg_', char(datetime)];
% FinalFolderName=strrep(SaveFolderName, ':', '.');
% mkdir(FinalFolderName) % Automatically generate a folder for saving
% 
% % Build a Global FALCON model for optimisation
% [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
% 
% % Define optimisation options
% estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
% estim.SSthresh=0.001;
% 
% estim.Lambda=2^-6;
% estim.Reg='L1/2';
% estim.Reg='none' 
% fxt_all=[];
% Costs_all=[];
% for cc=1:optRound
%     tic
%     k=FalconIC(estim,'normal'); %initial conditions
%     [xval,fval,AIC,MSE, BIC,Np]=FalconObjFun(estim,k); %objective function
%     Fitting_Time=toc;
% 
%     fxt_all=[fxt_all;[fval xval Fitting_Time]];
%     Costs_all=[Costs_all;[AIC, MSE, BIC,Np]];
%     beep; pause(0.5); beep;
% end
% 
% [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
% [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
% 
% save('L-Noreg')
% 
% %%
% PowerLambda=-10:0.5:10;
% ListLambda=[0,10.^PowerLambda];
% M=[];
% 
% for L=ListLambda
%     % Create a save folder
%     SaveFolderName=['Results_', ExpName, '_L12reg_', char(datetime)];
%     FinalFolderName=strrep(SaveFolderName, ':', '.');
%     mkdir(FinalFolderName) % Automatically generate a folder for saving
% 
%     % Build a Global FALCON model for optimisation
%     [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
% 
%     % Define optimisation options
%     estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
%     estim.SSthresh=0.001;
% 
%     estim.Lambda=L;
%     estim.Reg='L1/2';
%     fxt_all=[];
%     Costs_all=[];
%     for cc=1:optRound
%         tic
%         k=FalconIC(estim,'normal'); %initial conditions
%         [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
%         Fitting_Time=toc;
% 
%         fxt_all=[fxt_all;[fval xval Fitting_Time]];
%         Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
%         beep; pause(0.5); beep;
%     end
% 
%     [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
%     [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
%     idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
%     M=[M;fxt_all(idx,:),Costs_all(idx,:)];
% end
% 
% figure,
% subplot(3,1,1)
% plot(M(:,2:end-5)-M(1,2:end-5))
% set(gca,'XTick',1:length(ListLambda))
% set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% set(gca,'xticklabelrotation',90)
% subplot(3,1,2)
% plot(M(:,[2:5,34:37]))
% set(gca,'XTick',1:length(ListLambda))
% set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% set(gca,'xticklabelrotation',90)
% subplot(3,1,3)
% bar(M(:,43))
% set(gca,'XTick',1:length(ListLambda))
% set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% set(gca,'xticklabelrotation',90)
% legend('AIC')
% figure, bar(M(:,end))
% set(gca,'XTick',1:length(ListLambda))
% set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% set(gca,'xticklabelrotation',90)
% figure,
% plot(M(:,[]))
% 
% 
% 
% save('L12Large')
% %%
% PowerLambda=-20:1:10;
% ListLambda=[0,2.^PowerLambda];
% M=[];
% 
% for L=ListLambda
%     % Create a save folder
%     SaveFolderName=['Results_', ExpName, '_L12reg_', char(datetime)];
%     FinalFolderName=strrep(SaveFolderName, ':', '.');
%     mkdir(FinalFolderName) % Automatically generate a folder for saving
% 
%     % Build a Global FALCON model for optimisation
%     [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
% 
%     % Define optimisation options
%     estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
%     estim.SSthresh=0.001;
% 
%     estim.Lambda=L;
%     estim.Reg='L1/2';
%     fxt_all=[];
%     Costs_all=[];
%     for cc=1:optRound
%         tic
%         k=FalconIC(estim,'normal'); %initial conditions
%         [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
%         Fitting_Time=toc;
% 
%         fxt_all=[fxt_all;[fval xval Fitting_Time]];
%         Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
%         beep; pause(0.5); beep;
%     end
% 
%     [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
%     [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
%     idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
%     M=[M;fxt_all(idx,:),Costs_all(idx,:)];
% end
% 
% % figure,
% % subplot(3,1,1)
% % plot(M(:,2:end-5)-M(1,2:end-5))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % subplot(3,1,2)
% % plot(M(:,[2:5,34:37]))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % subplot(3,1,3)
% % bar(M(:,43))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % legend('AIC')
% % figure, bar(M(:,end))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % figure,
% % plot(M(:,[]))
% 
% save('L12Narrow')
% 
% 
% %% Single
% PowerLambda=-20:1:10;
% ListLambda=[0,2.^PowerLambda];
% M=[];
% 
% for L=ListLambda
%     % Create a save folder
%     SaveFolderName=['Results_', ExpName, '_L12reg_', char(datetime)];
%     FinalFolderName=strrep(SaveFolderName, ':', '.');
%     mkdir(FinalFolderName) % Automatically generate a folder for saving
% 
%     % Build a Global FALCON model for optimisation
%     [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
%     estim=FalconMakeModel(InputFile,char(MeasFileList(1)),0.5)
%     % Define optimisation options
%     estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
%     estim.SSthresh=0.001;
% 
%     estim.Lambda=L;
%     estim.Reg='L1/2';
%     fxt_all=[];
%     Costs_all=[];
%     for cc=1:optRound
%         tic
%         k=FalconIC(estim,'normal'); %initial conditions
%         [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
%         Fitting_Time=toc;
% 
%         fxt_all=[fxt_all;[fval xval Fitting_Time]];
%         Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
%         beep; pause(0.5); beep;
%     end
% 
%     [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
%     [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
%     idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
%     M=[M;fxt_all(idx,:),Costs_all(idx,:)];
% end
% 
% % save('L12_singleNarrow')
% 
% % figure,
% % subplot(4,1,1)
% % plot(M(:,2:end-5)-M(1,2:end-5))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % title('All parameters')
% % subplot(4,1,2)
% % plot(M(:,[2,10]))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % axis([0 length(ListLambda) 0 0.1])
% % title('k1/k1b')
% % subplot(4,1,3)
% % plot(M(:,[2,10]))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % title('k1/k1b')
% % axis([0 length(ListLambda) 0.9 1])
% % subplot(4,1,4)
% % bar(M(:,13))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % legend('AIC')
% % figure, bar(M(:,end))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % figure,
% % plot(M(:,[]))
% 
% save('L12_singleNarrow')
% %%
% %%% 
% PowerLambda=-10:0.5:20;
% ListLambda=[0,10.^PowerLambda];
% M=[];
% 
% for L=ListLambda
%     % Create a save folder
%     SaveFolderName=['Results_', ExpName, '_L1Groupreg_', char(datetime)];
%     FinalFolderName=strrep(SaveFolderName, ':', '.');
%     mkdir(FinalFolderName) % Automatically generate a folder for saving
% 
%     % Build a Global FALCON model for optimisation
%     [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
% 
%     % Define optimisation options
%     estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
%     estim.SSthresh=0.001;
% 
%     estim.Lambda=L;
%     estim.RegMatrix=(reshape(1:length(estim.param_vector),length(ContextsList), length(estim.param_vector)/length(ContextsList)))';
%     estim.Reg='L1Groups';
%     fxt_all=[];
%     Costs_all=[];
%     for cc=1:optRound
%         tic
%         k=FalconIC(estim,'normal'); %initial conditions
%         [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
%         Fitting_Time=toc;
% 
%         fxt_all=[fxt_all;[fval xval Fitting_Time]];
%         Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
%         beep; pause(0.5); beep;
%     end
% 
%     [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
%     [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
%     idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
%     M=[M;fxt_all(idx,:),Costs_all(idx,:)];
% end
% 
% % figure,
% % subplot(3,1,1)
% % plot(M(:,2:end-5)-M(1,2:end-5))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['   0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % subplot(3,1,2)
% % 
% % plot(M(:,[2:5,34:37]))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['   0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % subplot(3,1,3)
% % bar(M(:,43))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['   0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % legend('AIC')
% % figure, bar(M(:,end))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['   0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % figure,
% % plot(M(:,44))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['   0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% 
% save('LGroupsLarge')
% 
% 
% %%
% %%% 
% PowerLambda=-25:10;
% ListLambda=[0,2.^PowerLambda];
% M=[];
% 
% for L=ListLambda
%     % Create a save folder
%     SaveFolderName=['Results_', ExpName, '_L1Groupreg_', char(datetime)];
%     FinalFolderName=strrep(SaveFolderName, ':', '.');
%     mkdir(FinalFolderName) % Automatically generate a folder for saving
% 
%     % Build a Global FALCON model for optimisation
%     [estim] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound);
% 
%     % Define optimisation options
%     estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter); % Default setting
%     estim.SSthresh=0.001;
% 
%     estim.Lambda=L;
%     estim.RegMatrix=(reshape(1:length(estim.param_vector),length(ContextsList), length(estim.param_vector)/length(ContextsList)))';
%     estim.Reg='L1Groups';
%     fxt_all=[];
%     Costs_all=[];
%     for cc=1:optRound
%         tic
%         k=FalconIC(estim,'normal'); %initial conditions
%         [xval,fval,AIC,MSE,BIC,Np]=FalconObjFun(estim,k); %objective function
%         Fitting_Time=toc;
% 
%         fxt_all=[fxt_all;[fval xval Fitting_Time]];
%         Costs_all=[Costs_all;[AIC, MSE, BIC, Np]];
%         beep; pause(0.5); beep;
%     end
% 
%     [bestx,meanx,stdx]=FalconResults(estim,fxt_all,estim.param_vector,FinalFolderName);
%     [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,[PlotFitSummary PlotFitIndividual PlotHeatmapCost PlotStateSummary PlotStateEvolution],FinalFolderName);
%     idx=find(ismember(min(fxt_all(:,1)),fxt_all(:,1)));
%     M=[M;fxt_all(idx,:),Costs_all(idx,:)];
% end
% 
% % figure,
% % subplot(4,1,1)
% % plot(M(:,2:end-5)-M(1,2:end-5))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % title(' All Parameters')
% % 
% % subplot(4,1,2)
% % P=(M(:,2:end-5))';
% % P2=(std(reshape(P(:),4,length(P(:))/length(ContextsList))))';
% % S=reshape(P2, length(P2)/size(P,2), size(P,2));
% % Collapsed=sum(S<0.01);
% % Nparams=(max(Collapsed)-Collapsed).*length(ContextsList) + Collapsed;
% % bar(Nparams)
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % title('#Parameters')
% % 
% % subplot(4,1,3)
% % plot(M(:,6:9));
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % title('Example parameter k2')
% % legend({'Cell Line a','Cell Line b','Cell Line c','Cell Line d'})
% % 
% % subplot(4,1,4)
% % bar(M(:,43))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % title('AIC')
% % 
% % figure, bar(M(:,end))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % figure,
% % plot(M(:,44))
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)
% % 
% save('LGroupsNarrow')
% % 
% % 
% % figure, imagesc(P)
% % set(gca,'XTick',1:length(ListLambda))
% % set(gca,'xticklabels',['  0'; num2str(PowerLambda(:))])
% % set(gca,'xticklabelrotation',90)


%%
%%% 
PowerLambda1=-20:1:2; % L1/2
PowerLambda2=-20:1:2; % L1group
ListLambda1=[0,2.^PowerLambda1];
ListLambda2=[0,2.^PowerLambda2];

M=[];

for L1=ListLambda1
    for L2=ListLambda2
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

figure, 
AICs=M(:,end-3);
AICs=reshape(AICs,length(ListLambda1),length(ListLambda2))
imagesc(AICs)
title('AIC')
set(gca,'XTick',1:length(ListLambda1))
set(gca,'xticklabels',['  0'; num2str(PowerLambda2(:))])
set(gca,'xticklabelrotation',90)
set(gca,'YTick',1:length(ListLambda2))
set(gca,'yticklabels',['  0'; num2str(PowerLambda1(:))])
ylabel('L1 grouped')
xlabel('L1/2')

figure,
MSEs=M(:,end-2);
MSEs=reshape(MSEs,length(ListLambda1),length(ListLambda2))
imagesc(log(MSEs))
title('MSE')
ylabel('L1 grouped')
xlabel('L1/2')


figure, 
Nps=M(:,end);
Nps=reshape(Nps,length(ListLambda1),length(ListLambda2))
imagesc(Nps)
title('#Parameters')
ylabel('L1 grouped')
xlabel('L1/2')










