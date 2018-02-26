function [estim] = FalconKONodes(varargin)
% FalconKO creates new models for each knock-outed node by creating a dummy
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
bestx=varargin{2};
fxt_all=varargin{3};
MeasFile=varargin{4};
HLbound=varargin{5};
optRound_KO=varargin{6};
ToSave=0;
if nargin>6
    Folder=varargin{7};
    ToSave=1;
end

estim_orig = estim;
States_original=estim.state_names;
Param_original=estim.param_vector;
Interactions_original=estim.Interactions;
num_params=estim.param_vector;
bestcost = min(fxt_all(:,1));

%% AIC calculation
N = numel(estim.Output)-sum(sum(isnan(estim.Output)));
MSE= bestcost;
Nodes=estim.state_names;
Nodes(estim.Input_idx(1,:))=[];
pn=length(Nodes);
p= numel(Param_original);

AIC_complete = N*log(MSE) + 2*p; %AIC for base model

p_KD = zeros(1,pn);
param_vector=estim.param_vector;
SSthresh=estim.SSthresh;

cost_KD = zeros(1,pn);

%%% parameter perturbation and refitting
thisfig=figure; hold on
suptitle('Knock-out analysis (nodes)');

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
    estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',3000,'MaxIter',3000); % Default
    estim.SSthresh=1e-3;
    fval_all=[];
    
    for counterround=1:optRound_KO %computation for modified model
        k=FalconIC(estim); %initial conditions
        [xval,fval]=FalconObjFun(estim,k); %objective function
        fval_all=[fval_all; fval];
    end
    
    cost_KD(1,counter)=min(fval_all);
    
    %% reduced model (- parameter)
    
    N_r = numel(estim.Output)-sum(sum(isnan(estim.Output))); %number of datapoints
    p_r= numel(estim.param_vector); %number of parameters
    %remove parameters related to the knocked-out node
    Is=Interactions_original(strcmp(Interactions_original(:,2),thisNode),5); %fetch outgoing parameters
    Is=[Is;Interactions_original(strcmp(Interactions_original(:,4),thisNode),5)]; %fetch incoming parameters
    for cc=length(Is):-1:1 %remove numbers
        if str2num(char(Is(cc)))>=0
            Is(cc)=[];
        end
    end
    
    
    AIC_KD(counter) = N_r*log(cost_KD(counter)) + 2*(p_r-numel(unique(Is)));
    
    %%Plot AIC values
    
    AIC_merge=[AIC_complete,AIC_KD];
    AIC_merge_scaled = AIC_merge - AIC_complete   %% set ctrl model to 0
    set(0,'CurrentFigure',thisfig);
    figko=thisfig; hold on;
    
    h=bar(1,AIC_merge_scaled(1,1)); hold on
    for counter2 = 2:counter+1
        h=bar(counter2,AIC_merge_scaled(counter2)); hold on
        if AIC_merge_scaled(counter2) <= AIC_merge_scaled(1,1)
            set(h,'FaceColor','r');
        else%if AIC_merge(counter2) > AIC_complete
            set(h,'FaceColor','k');
        end
    end
    
    hline=refline([0 AIC_merge_scaled(1,1)]);
    hline.Color = 'r';
    
    set(gca,'XTick',[1:length(Nodes)+1])
    Xtitles=['Ctrl';Nodes'];
    set(gca, 'XTicklabel', Xtitles);
    %     title('AIC');
    xlabel('');
    set(gca, 'XTickLabelRotation', 45)
    ylabel('scaled Akaike Information Criterion (AIC)');
    hold off
    Min=min(AIC_merge_scaled(1:counter)); Max=max(AIC_merge_scaled(1:counter));
    if counter>1
        axis([0.5 counter+1.5 Min-0.1*abs(Max-Min) Max+0.1*abs(Max-Min)])
    end
    drawnow;
    if ToSave
        saveas(figko,[Folder,filesep,'Nodes Knock-Outs'],'fig')        
    end
    
end
close(wb);

if ToSave
    saveas(figko,[Folder,filesep,'Nodes Knock-Outs'],'tif')
    saveas(figko,[Folder,filesep,'Nodes Knock-Outs'],'fig')
    saveas(figko,[Folder,filesep,'Nodes Knock-Outs'],'jpg')
    saveas(figko,[Folder,filesep,'Nodes Knock-Outs'],'svg')
end

toc

estim = estim_orig;

estim.Results.KnockOutNodes.Parameters=Xtitles';
estim.Results.KnockOutNodes.AIC_values=AIC_merge;
estim.Results.KnockOutNodes.KO_effect=AIC_merge>=AIC_merge(1);
estim.Results.KnockOutNodes.Interpretation={'0 = no KO effect','1 = KO effect'};
estim.Results.KnockOutNodes.AIC_complete = AIC_complete;
estim.Results.KnockOutNodes.AIC_KD = AIC_KD;
delete('KDN_TempFile.txt')

end
