function [estim] = FalconKONodes(varargin)
% FalconKO creates new models for each knocked-out node by creating a dummy
% inhibiting node. Then calculates the fitness of the model and compares
% all KO models with AIC.
% [estim] = FalconKONodes(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO,FinalFolderName)
%
% :: Input values ::
% estim             complete model definition
% bestx             the vector of best optimised parameters
% fxt_all           all fitting costs, parameter values and running time during the different optimisations
% MeasFile          experimental data
% HLbound           qualitative threshold between high and low range of parameter values
% optRound_KO       number of optimization rounds
% FinalFolderName   name of the folder for saving results
%
% :: Output value(s) ::
% estim             updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

%fetching values from arguments
estim=varargin{1};
%%% TODO: remove argument 2 (redundant)
fxt_all=varargin{3};
HLbound=varargin{5};
optRound_KO=varargin{6};
ToSave=0;
if nargin>6
    Folder=varargin{7};
    ToSave=1;
end

estim_orig = estim;
Param_original=estim.param_vector;
Interactions_original=estim.Interactions;
MSE = min(fxt_all(:,1));

%% BIC calculation
N = numel(estim.Output)-sum(sum(isnan(estim.Output)));
Nodes=estim.state_names;
Nodes(estim.Input_idx(1,:))=[];
pn=length(Nodes);
p= numel(Param_original);

BIC_complete = N*log(MSE) + (log(N))*p; %BIC for base model

p_KD = zeros(1,pn);
PreviousOptions=estim.options;
SSthresh=estim.SSthresh;

cost_KD = zeros(1,pn);
cost_error = zeros(1,p);

%%% parameter perturbation and refitting
thisfig=figure; hold on
set(gca,'TickLabelInterpreter','none')
suptitle('Virtual Node KO');

wb = waitbar(0,'Please wait...');
for counter =  1:size(p_KD,2)
    waitbar(counter/pn,wb,sprintf('Running Nodes Knock-out Round %d out of %d ...',counter,pn))
    
    Interactions=Interactions_original;
    estim=estim_orig;
    thisNode=Nodes(counter);
    
    Interactions=[Interactions; ['ix', 'X_INHIB_X', '-|', thisNode, '1', 'N', 'D']];
    estim.state_names=[estim.state_names, 'X_INHIB_X'];
    
    estim.Input_idx=[estim.Input_idx,ones(size(estim.Input_idx,1),1).*length(estim.state_names)];
    estim.Input=[estim.Input,ones(size(estim.Input,1),1)];
    
    MeasFile=FalconData2File(estim);
    FalconInt2File(Interactions,'KDN_TempFile.txt')
    
    estim=FalconMakeModel('KDN_TempFile.txt',MeasFile,HLbound);
    estim.options = PreviousOptions;
    estim.SSthresh = SSthresh;
    fval_all=[];
    
    for counterround=1:optRound_KO %computation for modified model
        k=FalconIC(estim); %initial conditions
        [xval,fval]=FalconObjFun(estim,k); %objective function
        fval_all=[fval_all; fval];
    end
    
    cost_KD(1,counter)=min(fval_all);
    cost_error(counter) = std(fval_all);
    
    %% reduced model (- parameter)
    
    N_r = numel(estim.Output)-sum(sum(isnan(estim.Output))); %number of datapoints
    p_r= numel(estim.param_vector); %number of parameters
    %remove parameters related to the knocked-out node
    Is=Interactions_original(strcmp(Interactions_original(:,2),thisNode),5); %fetch outgoing parameters
    Is=[Is;Interactions_original(strcmp(Interactions_original(:,4),thisNode),5)]; %fetch incoming parameters
    for cc=length(Is):-1:1 %remove numbers
        if str2double(char(Is(cc)))>=0
            Is(cc)=[];
        end
    end
    
    
    BIC_KD(counter) = N_r*log(cost_KD(counter)) + (log(N_r))*(p_r-numel(unique(Is)));
    BIC_alt1=N_r*log(cost_KD(counter)+cost_error(counter)) + (log(N_r))*p_r;
    BIC_alt2=N_r*log(cost_KD(counter)-cost_error(counter)) + (log(N_r))*p_r;
    BIC_error(counter) = max(abs(BIC_KD(counter)-BIC_alt1),abs(BIC_KD(counter)-BIC_alt2));
    %%Plot BIC values
    
    BIC_merge=[BIC_complete,BIC_KD];
    BIC_Error_merge=[0,BIC_error];
    set(0,'CurrentFigure',thisfig);
    figko=thisfig; hold on;
    
    bar(BIC_complete, 'FaceColor', [0.95 0.95 0.95]); hold on
    for counter2 = 2:counter+1
        h=bar(counter2,BIC_merge(counter2)); hold on
        if BIC_merge(counter2) <= BIC_complete
            set(h,'FaceColor','g');
        else 
            set(h,'FaceColor','r');
        end
        errorbar(counter2,BIC_merge(counter2), BIC_Error_merge(counter2),'Color', 'k', 'LineWidth', 1); hold on
    end    
    
    plot([0 size(BIC_merge,2)+1],[BIC_complete BIC_complete],'-r'), hold on
    
    set(gca,'XTick',[1:length(Nodes)+1])
    Xtitles=['Full';Nodes'];
    set(gca, 'XTicklabel', Xtitles);
    xlabel('');
    set(gca, 'XTickLabelRotation', 45)
    ylabel('Bayesian Information Criterion (BIC)');
    hold off
    Min=min(BIC_merge(1:counter)-BIC_Error_merge(1:counter)); Max=max(BIC_merge(1:counter)+BIC_Error_merge(1:counter));
    axis([0.5 counter+1.5 Min-0.1*abs(Min) Max+0.1*abs(Max)])
    drawnow;
    if ToSave
        saveas(figko,[Folder,filesep,'KO_Nodes'],'fig')        
    end
    
end
close(wb);

if ToSave
    saveas(figko,[Folder,filesep,'KO_Nodes'],'tif')
    saveas(figko,[Folder,filesep,'KO_Nodes'],'fig')
    saveas(figko,[Folder,filesep,'KO_Nodes'],'jpg')
    saveas(figko,[Folder,filesep,'KO_Nodes'],'svg')
    saveas(figko,[Folder,filesep,'KO_Nodes'],'pdf')
end

toc

estim = estim_orig;

estim.Results.KnockOutNodes.Parameters=Xtitles';
estim.Results.KnockOutNodes.BIC_values=BIC_merge;
estim.Results.KnockOutNodes.KO_effect=BIC_merge>=BIC_merge(1);
estim.Results.KnockOutNodes.Interpretation={'0 = no KO effect','1 = KO effect'};
estim.Results.KnockOutNodes.BIC_complete = BIC_complete;
estim.Results.KnockOutNodes.BIC_KD = BIC_KD;
delete('KDN_TempFile.txt')

end
