function [Ks,Costs,estim] = FalconResample(varargin)
% FalconResample creates a new dataset based on the mean and standard deviation of the original dataset.
% Useful to assess the variability of the parameter estimates, in conjunction with FalconLPSA
% [Ks,Costs,estim]=FalconResample(estim, bestx, optRound_Sampling, New_Data_Sampling, CV_cutoff, FinalFolderName)
%
% :: Input values ::
% estim               complete model definition
% bestx               the vector of optimal parameters
% optRound_Sampling   the number of times the resampling of experimental replicates is performed
% New_Data_Sampling   the number of experimental replicates from the original data
% CV_cutoff           the coefficient of variation (CV) cutoff value to determine if the optimised parameters have a large distribution
% FinalFolderName     name of the folder for saving results
%
% :: Output values ::
% Ks                  optimised parameter values from the optmisation based on the resampled data
% Costs               optimised fitting costs from the optmisation based on the resampled data
% estim               updated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

estim=varargin{1};
bestx=varargin{2};
Rep=varargin{3};
N=varargin{4};

estim_orig=estim;

ToSave=0;

if nargin>4
    CV_cutoff=varargin{5};
end
if nargin>5
    Folder=varargin{6};
    ToSave=1;
end

Ks=[];
Costs=[];
Data=estim.Output;
SD=estim.SD;

h = waitbar(0,'Please wait...');
for exp=1:Rep
    new=[];
    waitbar(exp/Rep,h,sprintf('Running Resampling Round %d out of %d ...',exp,Rep))
    for n=1:N
        new(:,:,n)=normrnd(Data,SD);
    end
    
    estim.Output=mean(new,3);
    estim.SD=std(new,0,3);
    k=FalconIC(estim);
    [xval,fval]=FalconObjFun(estim,k);
    Ks=[Ks;xval];
    Costs=[Costs;fval];
end
close(h);

estim=estim_orig;

if nargin>4
    
    disp(['Cost for perturbated measurements: ', num2str(mean(Costs)), ' +- ', num2str(std(Costs))])
    h1= figure; hold on;
    errorbar(mean(Ks),std(Ks),'.b','LineWidth',2); hold on;
    ylim([0 1])
    plot(bestx,'*r','MarkerSize',12);
    title('Distribution of optimised parameters after resampling')
    hold on;
    set(gca,'Xtick',1:length(mean(Ks)),'XTickLabel',estim.param_vector,'XGrid','on');
    legend('RangeResampling','BestParam')
    
    if ToSave
        saveas(h1,[Folder, filesep, 'Resampling'],'tif');
        saveas(h1,[Folder, filesep, 'Resampling'],'fig');
        saveas(h1,[Folder, filesep, 'Resampling'],'jpg');
        saveas(h1,[Folder, filesep, 'Resampling'],'svg');
        
    end
    
    estim.Results.Resampling.Parameters=estim.param_vector';
    estim.Results.Resampling.OptimisedParameter=Ks;
    estim.Results.Resampling.OptimisedSD=std(Ks);
    estim.Results.Resampling.LargeSD=(std(Ks)./mean(Ks))*100>CV_cutoff;
    estim.Results.Resampling.Costs=Costs;
    
    Heading=cell(1,3);
    Heading(1,1)={'parameters'};
    Heading(1,2)={'mean'};
    Heading(1,3)={'S.D.'};
    
    Resampling=[estim.param_vector num2cell(mean(Ks,1)') num2cell(std(Ks,0,1)')];
    
    disp('Summary of resampling:')
    disp(' ')
    disp([Heading;Resampling])
    disp(' ')
        
end

if length(varargin)>5
    xlswrite([pwd filesep Folder filesep 'Summary_Resampling.xls'],[Heading;Resampling])
end
end