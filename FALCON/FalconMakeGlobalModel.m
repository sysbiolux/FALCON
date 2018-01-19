function [estim, varargout] = FalconMakeGlobalModel(InputFile,GlobalEdgesList,MeasFileList,ContextsList,HLbound)
% Allows for the expansion of logical networks from one seed network.
% InputFile is a interaction list (txt or xlsx). See 'FalconMakeModel'
% FixedEdgesList is an interaction list of only the edges that should
% remain fixed across contexts. If no edge should remain fixed just
% indicate one of the input edges or one that is fixed to 1 (cannot be
% empty)
% MeasFileList is a list (cell) of measuremnt files, see 'FalconMakeModel'
% ContextList is a cell with names of cell lines, time points, etc..
% HLBound, Forced, see 'FalconMakeModel'

estim=FalconMakeModel(InputFile,cell2mat(MeasFileList(1)),HLbound);


%%% Reading the fixed interactions file
FixedInteractions=[];
if ~isempty(GlobalEdgesList)
    Point=find(ismember(GlobalEdgesList,'.'),1,'last'); %finding the last point in the file name
    Ext=GlobalEdgesList(Point+1:end); %retrieving the extension
    if strcmp(Ext,'txt') %if text file
        fid=fopen(GlobalEdgesList,'r'); %open the file
        LineCounter=0; %zero the line counter
        while 1 %get out only when line is empty
            tline = fgetl(fid); %read a line
            if ~ischar(tline), break, end %break out of the loop if line is empty
            LineCounter=LineCounter+1; %count the lines
%             disp(tline) %display the line
            Input = regexp(tline,'\t','split'); %find the tabs in the line's text
            if length(Input)==5
                if ~exist('FixedInteractions')
                    FixedInteractions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))};
                else
                    FixedInteractions=[FixedInteractions; {['i' num2str(LineCounter)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))}];
                end
            elseif length(Input)==6
                if ~exist('FixedInteractions')
                    FixedInteractions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5)),cell2mat(Input(6))};
                else
                    FixedInteractions=[FixedInteractions; {['i' num2str(LineCounter)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5)),cell2mat(Input(6))}];
                end
            end
        end
        fclose(fid);
    elseif strcmp(Ext,'xls') || strcmp(Ext,'xlsx')
        [~,~,Other]=xlsread(GlobalEdgesList,1);
        LineCounter=2; %initialize the line counter
        while LineCounter<=size(Other,1) %get out only when line is empty
            Input = Other(LineCounter,:); %read a line
%             disp(Input) %display the line
            if length(Input)==5
                if ~exist('FixedInteractions')
                    FixedInteractions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5))};
                else
                    FixedInteractions=[FixedInteractions; {['i' num2str(LineCounter-1)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5))}];
                end
            elseif length(Input)==6
                if ~exist('FixedInteractions')
                    FixedInteractions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5)),cell2mat(Input(6))};
                else
                    FixedInteractions=[FixedInteractions; {['i' num2str(LineCounter-1)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5)),cell2mat(Input(6))}];
                end
            end
            LineCounter=LineCounter+1; %count the lines
        end
        
    end
end
FixParams=unique(FixedInteractions(:,5));

ListInputs=estim.state_names(unique(estim.Input_idx(:)));
Ncells=size(MeasFileList,2);

I=estim.Interactions;

Itot={};

for c=1:size(I,1) %for each interaction
    thisI=I(c,:);
    Itot=[Itot;repmat(thisI,Ncells,1)];
    
    if ~ismember(thisI(2),ListInputs)
        rep={};
        for r=1:Ncells
            rep=[rep; [char(thisI(2)),'-',char(ContextsList(r))]];
        end
        Itot(end-(Ncells-1):end,2)=rep;
    end
    
    if ~ismember(thisI(4),ListInputs)
        rep={};
        for r=1:Ncells
            rep=[rep; [char(thisI(4)),'-',char(ContextsList(r))]];
        end
        Itot(end-(Ncells-1):end,4)=rep;
    end
    
    if str2num(char(thisI(5)))>0
        
    elseif ~ismember(thisI(5),FixParams)
        rep={};
        for r=1:Ncells
            rep=[rep; [char(thisI(5)),'-',char(ContextsList(r))]];
        end
        Itot(end-(Ncells-1):end,5)=rep;
    end
end


for m=1:length(MeasFileList)
    thisMeasFile=char(MeasFileList(m));
    Point=find(ismember(thisMeasFile,'.'),1,'last'); %finding the last point in the file name
    Ext=thisMeasFile(Point+1:end); %retrieving the extension
    if strcmp(Ext,'txt') %if text file
        fid=fopen(thisMeasFile,'r');
        LineCounter=0;
        
        Input_vector=[];
        Input_index=[];
        Output_vector=[];
        Output_index=[];
        SD_vector=[];
        Annotation=[];
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            
            LineCounter=LineCounter+1;
            ReadIO = regexp(tline,'\t','split');
            Annotation_data=ReadIO(1);
            InputRaw=ReadIO(2);
            OutputRaw=ReadIO(3);
            ReadInput=strsplit(char(InputRaw),',');
            ReadOutput=strsplit(char(OutputRaw),',');
            
            if length(ReadIO)==3 % Shared index between output value and SD
                SDRaw=ReadIO(3);
                ReadSD=strsplit(char(SDRaw),',');
            end
            
            Annotation = [Annotation; Annotation_data];
            
            
            count_input=1;
            Input_idx_collect=[];
            Input_value_collect=[];
            for counter=1:(size(ReadInput,2)/2)
                idx_ReadInput=find(ismember(estim.state_names,ReadInput(count_input)));
                value_ReadInput=str2num(cell2mat(ReadInput(count_input+1)));
                Input_idx_collect=[Input_idx_collect idx_ReadInput];
                Input_value_collect=[Input_value_collect value_ReadInput];
                count_input=count_input+2;
            end
            
            Input_vector=[Input_vector; Input_value_collect];
            Input_index=[Input_index; Input_idx_collect];
            
            count_output=1;
            Output_idx_collect=[];
            Output_value_collect=[];
            for counter=1:(size(ReadOutput,2)/2)
                idx_ReadOutput=find(ismember(estim.state_names,ReadOutput(count_output)));
                value_ReadOutput=str2num(cell2mat(ReadOutput(count_output+1)));
                Output_idx_collect=[Output_idx_collect idx_ReadOutput];
                Output_value_collect=[Output_value_collect value_ReadOutput];
                count_output=count_output+2;
            end
            
            Output_vector=[Output_vector; Output_value_collect];
            Output_index=[Output_index; Output_idx_collect];
            
            if length(ReadIO)==4 % Shared index between output value and SD
                
                count_SD=1;
                SD_value_collect=[];
                for counter=1:(size(ReadSD,2)/2)
                    value_ReadSD=str2num(cell2mat(ReadSD(count_SD+1)));
                    SD_value_collect=[SD_value_collect value_ReadSD];
                    count_SD=count_SD+2;
                end
                
                SD_vector=[SD_vector; SD_value_collect];
                
            end
            
        end
        fclose(fid);
    elseif strcmp(Ext,'xls') || strcmp(Ext,'xlsx')
        [~,sheetnames] = xlsfinfo(thisMeasFile);
        [~,~,OtherIn]=xlsread(thisMeasFile,sheetnames{1});
        [~,~,OtherOut]=xlsread(thisMeasFile,sheetnames{2});
        [~,~,OtherErr]=xlsread(thisMeasFile,sheetnames{3});
        Input_index=[];
        Output_index=[];
        
        Input_vector=cell2mat(OtherIn(2:end,2:end));
        Annotation = (OtherIn(2:end,1));
        for jj=2:length(OtherIn(:,2))   %%%% modif Philippe OtherIn(:,1)
            Input_index_coll=[];
            for j=1:length(OtherIn(1,:))
                %                 if ~isnan(cell2mat(OtherIn(jj,j)))
                Input_index_coll=[Input_index_coll,find(ismember(estim.state_names,OtherIn(1,j)))];
                %                 else
                %                     Input_index_coll=[Input_index_coll,NaN]
                %                 end
            end
            Input_index=[Input_index;Input_index_coll];
        end
        OtherOut(strcmp(OtherOut,'NaN'))={NaN};
        OtherOut(strcmp(OtherOut,''))={NaN};
        Output_vector=cell2mat(OtherOut(2:end,:));
        for jj=2:length(OtherOut(:,1))
            Output_index_coll=[];
            for j=1:length(OtherOut(1,:))
                %                 if ~isnan(cell2mat(OtherOut(jj,j)))
                Output_index_coll=[Output_index_coll,find(ismember(estim.state_names,OtherOut(1,j)))];
                %                 else
                %                     Output_index_coll=[Output_index_coll,NaN];
                %                 end
            end
            Output_index=[Output_index;Output_index_coll];
        end
        OtherErr(strcmp(OtherErr,'NaN'))={NaN};
        OtherErr(strcmp(OtherErr,''))={NaN};
        if ~isempty(OtherErr)
            SD_vector=cell2mat(OtherErr(2:end,:));
        end
    end
    InputNames=estim.state_names(unique(Input_index,'stable'));
    InputValues=Input_vector;
    OutputNames=estim.state_names(unique(Output_index(~isnan(Output_index)),'stable'));
    for l=1:size(OutputNames,1)
        for c=1:size(OutputNames,2)
            OutputNames(l,c)=cellstr([char(OutputNames(l,c)),'-',char(ContextsList(m))]);
            
        end
    end
    
    InNamesSheet(:,:,m)=InputNames;
    OutNamesSheet(:,:,m)=OutputNames;
    InputsSheet(:,:,m)=InputValues;
    OutputsSheet(:,:,m)=Output_vector;
    SDSheet(:,:,m)=SD_vector;
    
    
end

Page1={};
Page1=[Page1;InNamesSheet(:,:,1)];
Page1=[Page1;num2cell(InputsSheet(:,:,1))];

Page2={};
Page3={};
for m=1:length(MeasFileList)
    temp=OutNamesSheet(:,:,m);
    temp=[temp;num2cell(OutputsSheet(:,:,m))];
    Page2=[Page2,temp];
    temp=OutNamesSheet(:,:,m);
    temp=[temp;num2cell(SDSheet(:,:,m))];
    Page3=[Page3,temp];
end
stamp=mat2str((floor(now*100000000)));
tempdir = tempname;
mkdir(tempdir);
tempfile=[tempdir filesep 'Results_' stamp '_.xls'];
varargout{1}=stamp;
varargout{2}=tempfile;

Page1nan=[zeros(1,size(Page1,2));cell2mat(cellfun(@isnan,Page1(2:end,:),'UniformOutput',0))];
Page1(Page1nan>0) ={'NaN'};

Page2nan=[zeros(1,size(Page2,2));cell2mat(cellfun(@isnan,Page2(2:end,:),'UniformOutput',0))];
Page2(Page2nan>0) ={'NaN'};

Page3nan=[zeros(1,size(Page3,2));cell2mat(cellfun(@isnan,Page3(2:end,:),'UniformOutput',0))];
Page3(Page3nan>0) ={'NaN'};

%Check for excel to determine which IO to use.
useexcel = isExcelPresent();

if useexcel
%     xlswrite(tempfile,Page1,1)
%     xlswrite(tempfile,Page2,2)
%     xlswrite(tempfile,Page3,3)
    try
        xlswrite(tempfile, ['Annotation'; Annotation] , 1, 'A1');
        xlswrite(tempfile, Page1 , 1, 'B1');
    %     xlswrite(tempfile, estim.Input, 1, 'B2');
        xlswrite(tempfile, Page2, 2, 'A1');
    %     xlswrite(tempfile, estim.Output, 2, 'A2');
        xlswrite(tempfile, Page3, 3, 'A1');
    %     xlswrite(tempfile, estim.SD, 3, 'A2');
    catch
        xlswrite(tempfile, ['Annotation'; Annotation] , 1, 'A1');
        xlswrite(tempfile, Page1 , 1, 'B1');
    %     xlswrite(tempfile, estim.Input, 1, 'B2');
        xlswrite(tempfile, Page2, 2, 'A1');
    %     xlswrite(tempfile, estim.Output, 2, 'A2');
        xlswrite(tempfile, Page3, 3, 'A1');
    %     xlswrite(tempfile, estim.SD, 3, 'A2');    
    end
else
    setupxlwrite();
%     xlswrite(tempfile,Page1,1)
%     xlswrite(tempfile,Page2,2)
%     xlswrite(tempfile,Page3,3)
    xlwrite(tempfile, ['Annotation'; Annotation] , 1, 'A1');
    xlwrite(tempfile, Page1 , 1, 'B1');
%     xlswrite(tempfile, estim.Input, 1, 'B2');
    xlwrite(tempfile, Page2, 2, 'A1');
%     xlswrite(tempfile, estim.Output, 2, 'A2');
    xlwrite(tempfile, Page3, 3, 'A1');
%     xlswrite(tempfile, estim.SD, 3, 'A2');
end



clearvars estim

GlobalFile=['GlobalInputFile_' stamp '.txt'];
FalconInt2File(Itot,GlobalFile);

estim=FalconMakeModel(GlobalFile,tempfile,HLbound);
disp(['There are ', num2str(size(MeasFileList,2)), ' different contexts'])

end
