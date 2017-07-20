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

estim=varargin{1};
fxt_all=varargin{2};
param_vector=varargin{3};
if length(varargin)>3
    FinalFolderName=varargin{4};
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
Heading_1stCol={'SSE';'Time(s)'};
FitCost_Time=[Heading_1stCol num2cell([bestfxt(1) bestfxt(end)]') num2cell([mean(fxt_all(:,1),1) mean(fxt_all(:,end),1)]') num2cell([std(fxt_all(:,1),0,1) std(fxt_all(:,end),0,1)]')];

disp('Summarized fitting cost and time:')
disp(' ')
disp([Heading;FitCost_Time])
disp(' ')

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
    useexcel = isExcelPresent();
    if useexcel        
        xlswrite([pwd filesep FinalFolderName filesep 'Summary_Optimised_Parameters.xls'],[Heading;param_name_best_mean_std])
    else
        tab = table(param_vector,bestx',meanx',stdx','VariableNames',{'Parameter','Best','Average','Std'});
        writetable(tab,[FinalFolderName, filesep, 'Summary_Optimised_Parameters.csv'],'Delimiter',',');
    end
end

%writing log file:
if exist(FinalFolderName)
    Start=FinalFolderName(9:end);
else
    Start='?';
end

Now=datestr(now);
fid = fopen([FinalFolderName, filesep, [FinalFolderName, '.log']],'at');
fprintf(fid, 'FALCON log file \n');
fprintf(fid, 'Start timestamp: %s \n', Start);
fprintf(fid, 'Finish timestamp: %s \n', Now);
fprintf(fid, 'Final MSE: %d \n', min(fxt_all(:,1)));
fprintf(fid, 'Network: \n');
[L,C]=size(estim.Interactions);
for c=1:L
    fprintf(fid, '%s %s %s %s %s %s %s \n', estim.Interactions{c,:});
end

fprintf(fid, 'Dataset: \n');
[L,C]=size(estim.Input);
FormatIn=repmat('%d ', 1, C);
fprintf(fid, 'Inputs: \n');
for c=1:L
    fprintf(fid, [FormatIn, '\n'], estim.Input(c,:));
end

[L,C]=size(estim.Input_idx);
FormatIn=repmat('%d ', 1, C);
fprintf(fid, 'Inputs indices: \n');
for c=1:L
    fprintf(fid, [FormatIn, '\n'], estim.Input_idx(c,:));
end

[L,C]=size(estim.Output);
FormatIn=repmat('%d ', 1, C);
fprintf(fid, 'Outputs: \n');
for c=1:L
    fprintf(fid, [FormatIn, '\n'], estim.Output(c,:));
end

[L,C]=size(estim.Output_idx);
FormatIn=repmat('%d ', 1, C);
fprintf(fid, 'Outputs indices: \n');
for c=1:L
    fprintf(fid, [FormatIn, '\n'], estim.Output_idx(c,:));
end

disp('==============================================')
disp(' ')

end