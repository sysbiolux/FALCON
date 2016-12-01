function [bestx, meanx, stdx] = FalconResults(varargin)
% FalconResults summarises the optimisation outputs and displays the their best and average values.
% [bestx,meanx,stdx]=FalconResults(fxt_all,param_vector,FinalFolderName)
%
% :: Input values ::
% fxt_all           all fitting costs, parameters and time during optimisations
% param_vector      parameter names (estim.param_vector)
% FinalFolderName   name of the folder for saving results
%
% :: Output values ::
% bestx             the vector of best optimised parameter values
% meanx             mean values of all optimised parameters
% stdx              standard deviations of all optimised parameters
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

fxt_all=varargin{1};
param_vector=varargin{2};
if length(varargin)>2
    FinalFolderName=varargin{3};
end

% Extract best parameter set and calculate means & SDs
bestx_index= find(ismember(fxt_all(:,1),min(fxt_all(:,1))));
bestfxt=fxt_all(bestx_index(1),:);
bestx=bestfxt(2:end-1);
meanx=mean(fxt_all(:,2:end-1),1);
stdx=std(fxt_all(:,2:end-1),0,1);

% Display results in table format

% First table: Fitting Cost and Time

Heading=cell(1,4);
Heading(1,1)={'-------'};
Heading(1,2)={'best'};
Heading(1,3)={'mean'};
Heading(1,4)={'S.D.'};
Heading_1stCol={'FitCost';'Time(s)'};
FitCost_Time=[Heading_1stCol num2cell([bestfxt(1) bestfxt(end)]') num2cell([mean(fxt_all(:,1),1) mean(fxt_all(:,end),1)]') num2cell([std(fxt_all(:,1),0,1) std(fxt_all(:,end),0,1)]')];

disp('Summarized fitting cost and time:')
disp(' ')
disp([Heading;FitCost_Time])
disp(' ')
if length(varargin)>2
end
% Second Table: Parameter values

% Convert parameter names from symbolic to string
param_string={};
for counter=1:size(param_vector,1)
    param_string(counter)={char(param_vector(counter))};
end
s=size(param_vector,1);

Heading=cell(1,4);
Heading(1,1)={'parameters'};
Heading(1,2)={'best'};
Heading(1,3)={'mean'};
Heading(1,4)={'S.D.'};

param_name_best_mean_std=[param_string' num2cell(bestfxt(2:end-1)') num2cell(meanx') num2cell(stdx')];

disp('Summarized identified parameters:')
disp(' ')
disp([Heading;param_name_best_mean_std])

if length(varargin)>2
    xlswrite([pwd filesep FinalFolderName filesep 'Summary_Optimised_Parameters.xls'],[Heading;param_name_best_mean_std]);
end

disp('==============================================')
disp(' ')

end