function [estim] = FalconKO(varargin)
% FalconKO creates new models for each knock-outed parameter by setting parameter value to 0 and re-optimise
% [estim] = FalconKO(estim, bestx, fxt_all, MeasFile, HLbound, optRound_KO,FinalFolderName)
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

Param_original=estim.param_vector;
Interactions_original=estim.Interactions;
num_params=estim.param_vector;
bestcost = min(fxt_all(:,1));

%% AIC calculation
N = numel(estim.Output);
SSE= bestcost;
p= numel(Param_original);

AIC_complete = N*log(SSE/N) + 2*p;

p_KD = zeros(1,p);
param_vector=estim.param_vector;
SSthresh=estim.SSthresh;


%%% Pick index and evaluate
cost_KD = zeros(1,p);

%%% parameter perturbation and refitting
thisfig=figure; hold on
suptitle('Knock-out analysis');

wb = waitbar(0,'Please wait...');
for counter =  1:size(p_KD,2)
    waitbar(counter/p,wb,sprintf('Running Knock-out Round %d out of %d ...',counter,p))
        
    %     for counter2 =  1:size(p_KD,1)
    Interactions=Interactions_original;
    replace_idx=find(ismember(Interactions_original(:,5),Param_original(counter)));
    disp('Removing interaction...')
    disp(Interactions_original(replace_idx,:))
    for counter2=1:length(replace_idx)
        Interactions{replace_idx(counter2),5}=char(num2str(p_KD(1, counter)));
    end
    
    FalconInt2File(Interactions,'KD_TempFile.txt')
    estim=FalconMakeModel('KD_TempFile.txt',MeasFile,HLbound,1);
    estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',3000,'MaxIter',1000); % Default
    estim.SSthresh=1e-3;
    fval_all=[];

    for counterround=1:optRound_KO
        k=FalconIC(estim); %initial conditions
        [xval,fval]=FalconObjFun(estim,k); %objective function
        fval_all=[fval_all; fval];
    end
    
    cost_KD(1,counter)=min(fval_all);
    
    %% reduced model (-1 parameter)
        
    N_r = numel(estim.Output);
    p_r= (numel(estim.param_vector));
    
    AIC_KD(counter) = N_r*log(cost_KD(counter)/N_r) + 2*p_r;
    
    %%Plot AIC values
    
    AIC_merge=[AIC_complete,AIC_KD];
    figko=figure(thisfig); hold on;
    
    h=bar(1,AIC_complete); hold on
    for counter2 = 2:counter+1
        h=bar(counter2,AIC_merge(counter2)); hold on
        if AIC_merge(counter2) <= AIC_complete
            set(h,'FaceColor','g');
        else%if AIC_merge(counter2) > AIC_complete
            set(h,'FaceColor','k');
        end
    end
    
    hline=refline([0 AIC_complete]);
    hline.Color = 'r';
    
    set(gca,'XTick',[1:length(num_params)+1])
    Xtitles=['Ctrl';Param_original];
    set(gca, 'XTicklabel', Xtitles);
    %     title('AIC');
    xlabel('');
    ylabel('Akaike Information Criterion (AIC)');
    hold off
    drawnow;
    if ToSave
        saveas(figko,[Folder,filesep,'Knock-Outs'],'tif')
        saveas(figko,[Folder,filesep,'Knock-Outs'],'fig')
        saveas(figko,[Folder,filesep,'Knock-Outs'],'jpg')

    end
    
end
close(wb);

toc

estim = estim_orig;

estim.Results.KnockOut.Parameters=Xtitles';
estim.Results.KnockOut.AIC_values=AIC_merge;
estim.Results.KnockOut.KO_effect=AIC_merge>=AIC_merge(1);


delete('KD_TempFile.txt')

end