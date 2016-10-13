function [cost_SA, estim]=FalconLPSA(varargin)
% FalconLPSA performs a local parameters sensitivity analysis (LPSA) on an optimized model
% [cost_SA, estim]=FalconLPSA(estim, bestx, MeasFile, HLbound, optRound_LPSA, increment, option, Parallisation, FinalFolderName);
%
% :: Input values ::
% estim             complete model definition
% bestx             the vector of best optimised parameters
% MeasFile          experimental data
% HLbound           qualitative threshold between high and low range of parameter values
% optRound_LPSA     number of optimization rounds for local parameter sensitivity analysis
% increment         number of parameter values' perturbation(s) spread from the optimised valued
% option            'complete' = pertubation over entire parameter space
%                   'fast = perturbation until reaching the assigned cut-off
%                   the assigned cut-off comes from mean+SD of fitting costs from FalconResampling
%
% :: Output values ::
% cost_SA           fitting cost for each perturbation
% estim             updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

estim=varargin{1};
bestx=varargin{2};
MeasFile=varargin{3};
HLbound=varargin{4};
optRound_LPSA=varargin{5};
LPSA_Increments=varargin{6};
ToSave=0;
if nargin>6
    IsFast=varargin{7};
    if nargin>7
        Parallelisation=varargin{8};
        if nargin>8
            Folder=varargin{9};
            ToSave=1;
        end
    end
end

estim_orig=estim;

%%%1) Resampling to choose cut-off
% if strcmp(IsFast,'fast')
[~,Costs]=FalconResample(estim,bestx, 3,10);
CutOff=mean(Costs)+std(Costs);
disp(['Cost for perturbated measurements: ', num2str(mean(Costs)), ' +- ', num2str(std(Costs))])
% end
SSthresh=estim.SSthresh;

p = bestx;
p_SA = zeros(2*LPSA_Increments+1,length(p));
p_SA(LPSA_Increments+1,:) = p;

for counter = 1:length(p)
    Min=estim.LB(counter);
    Max=estim.UB(counter);
    
    if ~isempty(estim.A)
        if sum(estim.A(:,counter))>0
            IdxC=find(estim.A(:,counter));
            Max=min(estim.b(IdxC),Max);
        end
    end
    
    for counter2 = 1:(LPSA_Increments-1)
        p_pert_values_pos =p(counter)+ ((Max-p(counter))/LPSA_Increments*counter2);
        p_pert_values_neg =p(counter)-((p(counter)-Min)/LPSA_Increments*counter2);
        p_SA(LPSA_Increments+1+counter2,counter) = p_pert_values_pos;
        p_SA(LPSA_Increments+1-counter2,counter) = p_pert_values_neg;
    end
end

%%% Pick index and evaluate
cost_SA = NaN(size(p_SA));
params_SA = NaN([size(p_SA),length(p)]);

%% Resimulation with perturbed parameter values
estim=estim_orig;
thisfig=figure; hold on
Param_original=estim.param_vector;
Interactions_original=estim.Interactions;
MidPoint=ceil(size(p_SA,1)/2);

num_plots=length(Param_original);
NLines=ceil(sqrt(num_plots));
NCols=ceil(num_plots/NLines);

Ident_All = [];

h = waitbar(0,'Please wait...');

for counter =  1:size(p_SA,2) %for each parameter
    waitbar(counter/size(p_SA,2),h,sprintf('Running LPSA Round %d out of %d ...',counter,size(p_SA,2)))
    
    OK=1;
    %1) for middle to 0
    for counter2 = MidPoint:-1:1
        if OK==1
            Interactions=Interactions_original;
            replace_idx=find(ismember(Interactions_original(:,5),Param_original(counter)));
            for counter3=1:length(replace_idx)
                Interactions{replace_idx(counter3),5}=char(num2str(p_SA(counter2, counter)));
            end
            FalconInt2File(Interactions,'LPSA_TempFile.txt');
            estim=FalconMakeModel('LPSA_TempFile.txt',MeasFile,HLbound,1);
            estim.Interactions
            estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',5000,'MaxIter',5000); % Default
            estim.SSthresh=SSthresh;
            
            toc_all=[];
            x_all=[];
            fval_all=[];
            
            if Parallelisation == 1
                
                parfor Round=1:optRound_LPSA %'parfor' will run the parallel computing toolbox
                    tic
                    k=FalconIC(estim,'normal'); %initial conditions
                    [xval,fval]=FalconObjFun(estim,k); %objective function
                    Fitting_Time=toc;
                    %collecting the run's outcome
                    toc_all=[toc_all; Fitting_Time];
                    x_all=[x_all; xval];
                    fval_all=[fval_all; fval];
                end
                
            else
                for Round=1:optRound_LPSA %'parfor' will run the parallel computing toolbox
                    tic
                    k=FalconIC(estim,'normal'); %initial conditions
                    [xval,fval]=FalconObjFun(estim,k); %objective function
                    Fitting_Time=toc;
                    %collecting the run's outcome
                    toc_all=[toc_all; Fitting_Time];
                    x_all=[x_all; xval];
                    fval_all=[fval_all; fval];
                end
            end
            
            fxt_all=[fval_all x_all toc_all];
            cost_SA(counter2,counter)=min(fval_all);
            params_SA(counter2,counter,find(ismember(Param_original,estim.param_vector)))=x_all(min(find(ismember(fval_all,min(fval_all)))),:)
            if strcmp(IsFast,'fast')
                if min(fval_all)>CutOff
                    OK=0;
                end
            end
        end
        
    end
    
    %2) from middle to 1
    OK=1;
    for counter2 =  MidPoint+1:size(p_SA,1)
        
        if OK==1
            Interactions=Interactions_original;
            replace_idx=find(ismember(Interactions_original(:,5),Param_original(counter)));
            for counter3=1:length(replace_idx)
                Interactions{replace_idx(counter3),5}=char(num2str(p_SA(counter2, counter)));
            end
            FalconInt2File(Interactions,'TempFile.txt');
            estim=FalconMakeModel('TempFile.txt',MeasFile,HLbound,1);
            estim.Interactions
            estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',5000,'MaxIter',5000); % Default
            estim.SSthresh=SSthresh;
            
            toc_all=[];
            x_all=[];
            fval_all=[];
            
            
            if Parallelisation == 1
                
                parfor Round=1:optRound_LPSA %'parfor' will run the parallel computing toolbox
                    tic
                    k=FalconIC(estim,'normal'); %initial conditions
                    [xval,fval]=FalconObjFun(estim,k); %objective function
                    Fitting_Time=toc;
                    %collecting the run's outcome
                    toc_all=[toc_all; Fitting_Time];
                    x_all=[x_all; xval];
                    fval_all=[fval_all; fval];
                end
                
            else
                for Round=1:optRound_LPSA %'parfor' will run the parallel computing toolbox
                    tic
                    k=FalconIC(estim,'normal'); %initial conditions
                    [xval,fval]=FalconObjFun(estim,k); %objective function
                    Fitting_Time=toc;
                    %collecting the run's outcome
                    toc_all=[toc_all; Fitting_Time];
                    x_all=[x_all; xval];
                    fval_all=[fval_all; fval];
                end
            end
            
            fxt_all=[fval_all x_all toc_all];
            
            cost_SA(counter2,counter)=min(fval_all);
            params_SA(counter2,counter,find(ismember(Param_original,estim.param_vector)))=x_all(min(find(ismember(fval_all,min(fval_all)))),:);
            
            if strcmp(IsFast,'fast')
                if min(fval_all)>CutOff
                    OK=0;
                end
            end
        end
        
    end
    
    
    % evaluate the identifiability of each parameter
    EvalSA_Temp = cost_SA(:,counter)<CutOff;
    EvalSA = ~EvalSA_Temp;
    Ident_Left  = sum(EvalSA(1:floor(length(EvalSA)/2)));
    Ident_Right = sum(EvalSA(ceil(length(EvalSA)/2)+1:length(EvalSA)));
    if sum([Ident_Left Ident_Right])==length(EvalSA)-1
        Identifiability = 1;
    elseif Ident_Left||Ident_Right > 0
        Identifiability = 2;
    elseif sum([Ident_Left Ident_Right])==0
        Identifiability = 3;
    end
    
    Ident_All(counter)=Identifiability;
    
    try
        % plot the results
        figure(thisfig), hold on,
        subplot(NLines,NCols,counter), hold on
        plot(p_SA(:,counter),cost_SA(:,counter), '-o','MarkerSize',5)
        hold on, plot(p_SA(LPSA_Increments+1,counter),cost_SA(LPSA_Increments+1,counter),'b*','MarkerSize',15)
        min_val=min(cost_SA(:,counter));
        %%% remove red box not to confuse user with optimum
        % min_index=find(min_val==cost_SA(:,counter));
        % hold on, plot(p_SA(min_index,counter),cost_SA(min_index,counter),'rs','MarkerSize',15)
        
        xlabel('Parameter range')
        ylabel('SSE')
        title(Param_original(counter))
        hline=refline([0 CutOff]);
        hline.Color = 'r';
        drawnow;
    catch
        warning('This figure is LPSA figure is not plotted')
    end
end

close(h);

try
    suptitle('Local parameter sensitivity analysis');
catch
end

if ToSave
    saveas(thisfig,[Folder, filesep, 'Identifiability'],'tif')
    saveas(thisfig,[Folder, filesep, 'Identifiability'],'fig')
    saveas(thisfig,[Folder, filesep, 'Identifiability'],'jpg')
    saveas(thisfig,[Folder, filesep, 'Identifiability'],'svg')
    
end


%% Parameter covariance

% figure, hold on
% Np=length(p);
% sub=1;
% for ppl=1:Np
%     for ppc=1:Np
%         subplot(Np,Np,sub)
%         thisfig2=plot(p_SA(:,ppl),params_SA(:,ppl,ppc));
%         axis([0 1 max(min(params_SA(:,ppl,ppc))-0.1,0) min(max(params_SA(:,ppl,ppc))+0.1,1)])
%         xlabel(Param_original(ppl))
%         ylabel(Param_original(ppc))
%         sub=sub+1;
%     end
% end
%
% suptitle('Covariance');
%
% if ToSave
%     saveas(thisfig2,[Folder, filesep, 'Covariance'],'tif')
%     saveas(thisfig2,[Folder, filesep, 'Covariance'],'fig')
%     saveas(thisfig2,[Folder, filesep, 'Covariance'],'jpg')
%     saveas(thisfig2,[Folder, filesep, 'Covariance'],'svg')

% end

estim=estim_orig;

estim.Results.LPSA.ParamNames=estim.param_vector';
estim.Results.LPSA.Identifiability=Ident_All;
estim.Results.LPSA.LPSA_Increments=LPSA_Increments;
estim.Results.LPSA.p_SA = p_SA;
estim.Results.LPSA.cost_SA = cost_SA;
estim.Results.LPSA.CutOff = CutOff;
estim.Results.LPSA.Interpretation={'1=Identifiable','2=Partially identifiable','3=Non-identifiable'};


Heading=cell(1,2);
Heading(1,1)={'parameters'};
Heading(1,2)={'Identifiable?'};


LPSA=[estim.param_vector num2cell(estim.Results.LPSA.Identifiability') ];

disp('Summary of local parameter sensistivity analysis:')
disp(' ')
disp([Heading;LPSA])
disp(' ')

if length(varargin)>8
    xlswrite([Folder filesep 'Summary_LPSA.xls'],[Heading;LPSA])
end



delete('LPSA_TempFile.txt')
delete('TempFile.txt')

end