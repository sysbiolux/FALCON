function [estim]=FalconMakeModel(InputFile,MeasFile,HLbound)
% FalconMakeModel creates the integrated optimization problem for Falcon from a list of interactions (InputFile)
% and a list of measured nodes in different conditions (MeasFile). The files can be .txt 
% (see 'example_model.txt' and 'example_meas.txt') or Excel (see 'PDGF.xlsx' and 'PDGF_meas.xlsx'). 
% Automatically assigns the weights of interactions for single-input node(s) and creates the constrains for the optimization problem.
% [estim]=FalconMakeModel(InputFile,MeasFile,HLbound,Forced)
% 
% :: Input values ::
% InputFile          model file (either txt or xls(x))
% MeasFile           experimental data file (either txt or xls(x))
% HLbound            qualitative threshold between high and low range of parameter values
% Forced             this variable has been removed 
%
% :: Output values::
% estim              integrated model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

clearvars Interactions

disp('...reading network file...')
%%% Reading the model file
Point=find(ismember(InputFile,'.'),1,'last'); %finding the last point in the file name
Ext=InputFile(Point+1:end); %retrieving the extension
if strcmp(Ext,'txt') || strcmp(Ext,'sif')%if text file
    fid=fopen(InputFile,'r'); %open the file
    LineCounter=0; %zero the line counter
    while 1 %get out only when line is empty
        tline = fgetl(fid); %read a line
        if ~ischar(tline), break, end %break out of the loop if line is empty
        LineCounter=LineCounter+1; %count the lines
%         disp(tline) %display the line
        Input = regexp(tline,'\t','split'); %find the tabs in the line's text
        if length(Input)==3 %if no param names, no gate, no high/low constrains
            Input{4}=[char(Input{1}),char(Input{2}),char(Input{3})]
            Input{5}='N'
            if ~exist('Interactions')
                Interactions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))};
            else
                Interactions=[Interactions; {['i' num2str(LineCounter)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))}];
            end
        elseif length(Input)==4 %if there are param names
            Input{5}='N'
            if ~exist('Interactions')
                Interactions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))};
            else
                Interactions=[Interactions; {['i' num2str(LineCounter)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))}];
            end
        elseif length(Input)==5 %if there are gates
            if ~exist('Interactions')
                Interactions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))};
            else
                Interactions=[Interactions; {['i' num2str(LineCounter)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5))}];
            end
        elseif length(Input)==6 %if there are high/low constrains
            if ~exist('Interactions')
                Interactions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5)),cell2mat(Input(6))};
            else
                Interactions=[Interactions; {['i' num2str(LineCounter)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),cell2mat(Input(4)),cell2mat(Input(5)),cell2mat(Input(6))}];
            end
        end        
    end
    fclose(fid);
elseif strcmp(Ext,'xls') || strcmp(Ext,'xlsx')
    [~,~,Other]=xlsread(InputFile,1);
    LineCounter=2; %initialize the line counter
    while LineCounter<=size(Other,1) %get out only when line is empty
        Input = Other(LineCounter,:); %read a line
%         disp(Input) %display the line
        if length(Input)==5 
            if ~exist('Interactions')
                Interactions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5))};
            else
                Interactions=[Interactions; {['i' num2str(LineCounter-1)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5))}];
            end
        elseif length(Input)==6
            if ~exist('Interactions')
                Interactions={'i1',cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5)),cell2mat(Input(6))};
            else
                Interactions=[Interactions; {['i' num2str(LineCounter-1)],cell2mat(Input(1)),cell2mat(Input(2)),cell2mat(Input(3)),num2str(cell2mat(Input(4))),cell2mat(Input(5)),cell2mat(Input(6))}];
            end
        end
        LineCounter=LineCounter+1; %count the lines
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QUALITY CONTROL %%%
disp('...checking network consistency...')


state_names=unique([Interactions(:,2);Interactions(:,4)]'); % extract state names

for counter=1:size(state_names,2) % for each state
    output_index=find(ismember(Interactions(:,4),state_names(counter))); % get the indices of the outputs of the current state
    if size(output_index,1)>1 % if there exist more than one output
        AND_Index = char(Interactions(output_index,6)) == 'A'; % get the AND indices
        OR_Index  = char(Interactions(output_index,6)) == 'O'; % get the OR indices
        NG_Index  = char(Interactions(output_index,6)) == 'N'; % get the NO GATE indices
        if ~prod(NG_Index) % if the product of all NO GATE indicies is NOT FALSE (i.e. there exist at least one interaction with AND/OR gate assigned)
            if sum(AND_Index) > 1 || sum(OR_Index) > 1 % check if the sum of all AND or OR indicies is more than 1 (i.e. there exist more than one interaction with AND/OR gate)
                if length(unique(Interactions(output_index(AND_Index),5)))==1 || length(unique(Interactions(output_index(OR_Index),5)))==1 % check if the parameters of those AND/OR interactions are the same
                    if length(unique(Interactions(output_index(AND_Index),3)))>1 || length(unique(Interactions(output_index(OR_Index),3)))>1 % if the same parameter is assigned, check if they have the same type of interaction
                        error(['Error: Please check the interactions lines ' num2str(output_index')]) % Error: the types of AND/OR interactions are differed with the same gate/same parameter
                    end                        
                else
                    error(['Error: Please check the interactions lines ' num2str(output_index')]) % Error: the AND/OR interactions are assigned with different parameters
                end
            else
                error(['Error: Please check the interactions lines ' num2str(output_index')]) % Error: there exist only a single interaction with AND/OR gate (need a pair of interactions for AND/OR gate)
            end
        else % If all interactions have NO GATE 
            if length(unique(Interactions(output_index,3))) > 1 && length(unique(Interactions(output_index,5))) == 1 && ~strcmp(unique(Interactions(output_index,5)),'1') % check if the interactions have different types of interaction AND if they have the same parameters AND if they are NOT fixed to 1
                error(['Error: Please check the interactions lines ' num2str(output_index')]) % Error: there exist NO GATE interactions with different types of interaction but same parameters which were not assigned to 1
            end            
        end   
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NETWORK EXPANSION %%%

DidChange=1;
while DidChange==1
    [DidChange, InteractionsNew]=FalconExpand(Interactions);
    Interactions=InteractionsNew;
    if DidChange, disp('...network has been expanded...'), end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WEIGHTS ASSIGNMENT %%%
%%% Fixing sum of probabilities to 1 %%%

IsBool=~strcmp(Interactions(:,6),'N');
IdxAct=strcmp(Interactions(:,3),'->');
IdxInh=strcmp(Interactions(:,3),'-|');
IdxNAct=(~IsBool.*IdxAct);
IdxNInh=(~IsBool.*IdxInh);
ReactInNAct=Interactions(IdxNAct==1,:);
ReactInNInh=Interactions(IdxNInh==1,:);
UOutInNAct=unique(ReactInNAct(:,4),'stable');
UOutInNInh=unique(ReactInNInh(:,4),'stable');
state_names=unique([Interactions(:,2);Interactions(:,4)]'); % extract state names
ma=zeros(size(state_names,2));
mi=zeros(size(state_names,2));

Interactions(IsBool,5)=repmat({'1'},sum(IsBool),1);

for i=1:length(UOutInNAct) % for each Output node in single interactions
    thisOut=UOutInNAct(i); % node ID
    isAlone=strcmp(thisOut,ReactInNAct(:,4)); % vector of identities in activating links
    if sum(isAlone)==1 % if this is the only activating interaction
        Idx=(strcmp(thisOut,Interactions(:,4))) & IdxAct; % look for it
        Interactions(Idx,5)=repmat({'1'},sum(Idx),1); % fix to 1
    end
    if sum(isAlone)>1 % if there are several interactions
        SumFixed=0;
        Idx=(strcmp(thisOut,Interactions(:,4))) & IdxAct;
        thisParams=Interactions(Idx,5); % activating parameters related to this node
        fParams=[];
        for cf=1:length(thisParams) % for each one
            if str2num(char(thisParams(cf)))>=0 % if number
                fParams=[fParams;1]; % vector of boolean for fixed
                SumFixed=SumFixed+str2num(char(thisParams(cf))); % adding to the sum of probabilities
            else
                fParams=[fParams;0]; % vector of boolean
            end
        end
        if sum(fParams)==length(fParams)-1 % if there is only one remaining non-fixed interaction
            % fix the parameter to the rest of probabilities
            Interactions(find(ismember(Interactions(:,5),thisParams(~fParams))),5)={num2str(1-SumFixed)};
        end
        if SumFixed==1
            ParamsToMod=thisParams(~fParams);
            for ii=1:length(ParamsToMod)
                thisMod=ParamsToMod(ii);
                Interactions(find(ismember(Interactions(:,5),thisMod)),5)={'0'};
            end
        end
    end
    
end


for i=1:length(UOutInNInh) % for each Output node in inhibitions
    thisOut=UOutInNInh(i); % node ID
    isAlone=strcmp(thisOut,ReactInNInh(:,4)); % vector of identities in activating links
    if sum(isAlone)>1 % if there are several interactions
        SumFixed=0;
        Idx=(strcmp(thisOut,Interactions(:,4))) & IdxInh;
        thisParams=Interactions(Idx,5); % inhibiting parameters related to this node
        fParams=[];
        for cf=1:length(thisParams) % for each one
            if str2num(char(thisParams(cf)))>=0 % if number
                fParams=[fParams;1]; % vector of boolean for fixed
                SumFixed=SumFixed+str2num(char(thisParams(cf))); % adding to the sum of probabilities
            else
                fParams=[fParams;0]; % vector of boolean
            end
        end
        if SumFixed==1
            Interactions(find(ismember(Interactions(:,5),thisParams(find(~fParams)))),5)={'0'};
        end
    end
end

% disp(Interactions)

UniqueParamNames=unique(Interactions(:,5),'stable')'; %parameter names
for u=length(UniqueParamNames):-1:1 %removing numeric values, starts from the end (otherwise suppresses the wrong line)
    if str2num(char(UniqueParamNames(u)))>=0
        UniqueParamNames(u)=[];
    end
end

% At this point the "Interactions" is in its final state.
FalconInt2File(Interactions,'Expanded.txt')

for i=1:size(Interactions,1) % for each interaction
    idxout=find(ismember(state_names,Interactions(i,4)));
    idxin=find(ismember(state_names,Interactions(i,2)));
    pv=str2num(char(Interactions(i,5)));
    % write the value in ma or mi, or 1 if free parameter
    if strcmp(Interactions(i,3),'->')
        if pv>=0
            ma(idxout,idxin)=pv;
        else
            ma(idxout,idxin)=1;
        end
    else
        if pv>=0
            mi(idxout,idxin)=pv;
        else
            mi(idxout,idxin)=1;
        end
    end
end

delete('Expanded.txt')

% Adding auto-activation for nodes without input
NoInput=find(sum(ma,2)==0);
for ai=1:length(NoInput)
    ma(NoInput(ai),NoInput(ai))=1;
end

% Infering constrains and making param index
Aeq=[];beq=[];A=[];b=[];

param_index=[];
cntb=0;
UniqueOutputList=unique(Interactions(:,4)); %list of outputs
FixBool=[];PN={};

for i=1:length(UniqueOutputList)%for each unique output
    thisOut=UniqueOutputList(i);
    Idx=find(ismember(Interactions(:,4),thisOut));
    GateType=unique(Interactions(Idx,6));
    if ~strcmp(GateType,'N') %if Bool gate
        isAct=strcmp(unique(Interactions(Idx,3)),'->');
        pGate=strcmp(GateType,'O')+strcmp(GateType,'A')*2; %type of boolean interaction (1=OR, 2=AND)
        pHL=0; %default HL bound
        if size(Interactions,2)==7 %only if 7 columns
            pHL=strcmp(Interactions(Idx,7),'L')*(-1)+strcmp(Interactions(Idx,7),'H'); % (-1=L, 1=H)
        end
        cntb=cntb+1;
        idxo=find(ismember(state_names,thisOut)); %get the index of that output
        idxi=find(ismember(state_names,Interactions(Idx,2)));
        param_index=[param_index; idxo, idxi(1), isAct, ~isAct, cntb, pGate, pHL(1)];
        param_index=[param_index; idxo, idxi(2), isAct, ~isAct, cntb, pGate, pHL(2)];
        FixBool=[FixBool,1,1]; PN=[PN,1,1];
        
    else % if only single interactions
        for ii=1:length(Idx) % for each interactions
            thisIdx=Idx(ii); % its index
            if str2num(char(Interactions(thisIdx,5)))>=0 % if number: not going to param_index (but is in ma and mi)
            else
                idxo=find(ismember(state_names,thisOut)); %get the index of that output
                idxi=find(ismember(state_names,Interactions(thisIdx,2))); %index of the input
                isAct=strcmp(unique(Interactions(thisIdx,3)),'->'); %type of interaction
                pHL=0; %default HL bound
                if size(Interactions,2)==7 %only if 7 columns
                    pHL=strcmp(Interactions(thisIdx,7),'L')*(-1)+strcmp(Interactions(thisIdx,7),'H'); % (-1=L, 1=H)
                end
                param_index=[param_index; idxo, idxi, isAct, ~isAct, 0, 0, pHL];
                FixBool=[FixBool,0]; PN=[PN,Interactions(thisIdx,5)];
                
            end
            
        end
        % get indices relating to that output
        IdxAct=find(ismember(strcmp(Interactions(:,3),'->').*strcmp(Interactions(:,4),thisOut),1));
        IdxInh=find(ismember(strcmp(Interactions(:,3),'-|').*strcmp(Interactions(:,4),thisOut),1));
        % fetch the names of the parameters, remove numeric values
        ParamAct=Interactions(IdxAct,5);
        SumAct=0;
        for kp=length(ParamAct):-1:1
            if str2num(char(ParamAct(kp)))>=0
                SumAct=SumAct+str2num(char(ParamAct(kp)));
                ParamAct(kp)=[];
            end
        end
        ParamInh=Interactions(IdxInh,5);
        SumInh=0;
        for kp=length(ParamInh):-1:1
            if str2num(char(ParamInh(kp)))>=0
                SumInh=SumInh+str2num(char(ParamInh(kp)));
                ParamInh(kp)=[];
            end
        end
        if numel(ParamAct)>1
            Aeq=[Aeq;zeros(1,length(UniqueParamNames))];
            idxadd=find(ismember(UniqueParamNames,ParamAct)); %get indexes of the additive inputs
            Aeq(end,idxadd)=1; %these parameters have to sum up to...
            beq=[beq; 1-SumAct]; %...remaining
        end
        if numel(ParamInh)>1
            A=[A;zeros(1,length(UniqueParamNames))];
            idxadd=find(ismember(UniqueParamNames,ParamInh)); %get indexes of the additive inputs
            A(end,idxadd)=1; %these parameters have to be smaller than...
            b=[b; 1-SumInh]; %...remaining
        end
        
    end
    
end

% derive additional High/Low constraints with bounds

LB=zeros(1,length(UniqueParamNames)); %create lower boundaries with default=0
UB=ones(1,length(UniqueParamNames)); %create upper boundaries with default=1

if size(Interactions,2)==7 %if interaction file contains H/L specifications
    for p=1:length(UniqueParamNames) %for each parameter
        idxp=find(ismember(Interactions(:,5),UniqueParamNames(p))); %get indexes of this parameter
        UL=unique(Interactions(idxp,7)); %get the bounds for that parameter (should be only one)
        
        if strcmp(UL,'L') %if there are low parameters
            UB(p)=HLbound; %set upper boundaries as HLbound
        end
        if strcmp(UL,'H') %if there are high parameters
            LB(p)=HLbound; %set lower boundaries as HLbounds
        end
    end
end

param_vector=UniqueParamNames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the info from the measurement file
disp('...reading datafile...')
if size(MeasFile,1)==1
    Point=find(ismember(MeasFile,'.'),1,'last'); %finding the last point in the file name

    Ext=MeasFile(Point+1:end); %retrieving the extension
    if strcmp(Ext,'txt') %if text file
        fid=fopen(MeasFile,'r');
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
            Annotation_data=string(ReadIO(1));
            InputRaw=ReadIO(2);
            OutputRaw=ReadIO(3);
            ReadInput=strsplit(char(InputRaw),',');
            ReadOutput=strsplit(char(OutputRaw),',');

            if length(ReadIO)==4 % Shared index between output value and SD
                SDRaw=ReadIO(4);
                ReadSD=strsplit(char(SDRaw),',');
            end

    %        Annotation= [];
           Annotation=[Annotation; Annotation_data];


            count_input=1;
            Input_idx_collect=[];
            Input_value_collect=[];
            for counter=1:(size(ReadInput,2)/2)
                idx_ReadInput=find(ismember(state_names,ReadInput(count_input)));
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
                idx_ReadOutput=find(ismember(state_names,ReadOutput(count_output)));
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

        [~,sheetnames] = xlsfinfo(MeasFile);

        [~,~,OtherIn]=xlsread(MeasFile,sheetnames{1});
        [~,~,OtherOut]=xlsread(MeasFile,sheetnames{2});
        [~,~,OtherErr]=xlsread(MeasFile,sheetnames{3});
        Input_index=[];
        Output_index=[];

        Input_vector=cell2mat(OtherIn(2:end,2:end)); 
        Annotation = cellfun(@num2str,(OtherIn(2:end,1)),'UniformOutput',0); 
        OtherIn

        for jj=2:length(OtherIn(:,2))
            Input_index_coll=[];
            for j=1:length(OtherIn(1,:))
                Input_index_coll=[Input_index_coll,find(ismember(state_names,OtherIn(1,j)))];
            end
            Input_index=[Input_index;Input_index_coll];
        end
        OtherOut(strcmp(OtherOut,'NaN'))={NaN};
        OtherOut(strcmp(OtherOut,''))={NaN};

        Output_vector=cell2mat(OtherOut(2:end,:));
        for jj=2:length(OtherOut(:,1))
            Output_index_coll=[];
            for j=1:length(OtherOut(1,:))
                Output_index_coll=[Output_index_coll,find(ismember(state_names,OtherOut(1,j)))];
            end
            Output_index=[Output_index;Output_index_coll];
        end
        OtherErr(strcmp(OtherErr,'NaN'))={NaN};
        OtherErr(strcmp(OtherErr,''))={NaN};
        if ~isempty(OtherErr)
            SD_vector=cell2mat(OtherErr(2:end,:));
        end
    end
    
else
    Input_index=[];
    Output_index=[];
    data=readtable(char(MeasFile(1)));
    Input_vector=table2array(data(:,2:end));
    Input_names=data.Properties.VariableNames(2:end);
    Annotation=table2array(data(:,1));
    
    data=readtable(char(MeasFile(2)));
    Output_vector=table2array(data);
    Output_names=data.Properties.VariableNames;
    
    if size(MeasFile,1)>2
        data=readtable(char(MeasFile(3)));
        SD_vector=table2array(data);
    else
        SD_vector=zeros(size(Output_vector));
    end
    
    for jj=1:length(Input_vector(:,1))
        Input_index_coll=[];
        for j=1:length(Input_vector(1,:))
            Input_index_coll=[Input_index_coll,find(ismember(state_names,Input_names(j)))];
        end
        Input_index=[Input_index;Input_index_coll];
    end
    
    for jj=1:length(Output_vector(:,1))
        Output_index_coll=[];
        for j=1:length(Output_vector(1,:))
            Output_index_coll=[Output_index_coll,find(ismember(state_names,Output_names(j)))];
        end
        Output_index=[Output_index;Output_index_coll];
    end
    
    
end


NewOutput_index=[];
Nexp=size(Output_index,1);
for oi=1:size(Output_index,2)
    thisId=Output_index(:,oi);
    thisId(isnan(thisId))=[];
    NewOutput_index=[NewOutput_index,unique(thisId)];
end
Output_index=repmat(NewOutput_index,Nexp,1);


kInv=[]; kInd=[];
for i=1:length(PN)
    thisParam=PN(i);
    if str2num(char(cell2mat(thisParam)))>=0
    else
        kInv=[kInv; find(ismember(param_vector,char(cell2mat(thisParam))))];
    end
end
kInd=kInv;

BoolIdx=[];
BoolLines=param_index(:,6)>0;
BoolOuts=unique(param_index(BoolLines,1));
BoolLines2=ismember(param_index(:,1),BoolOuts);
BoolIdx=param_index(:,5).*BoolLines2;
BoolIdx(~BoolIdx)=[];
BoolIdx=unique(BoolIdx,'stable');

pd=param_index(~FixBool,:);
IdxInAct=pd(pd(:,3)>0,2); %indices of the inputs in activating links
IdxOutAct=pd(pd(:,3)>0,1); %indices of the outputs in activating links
IdxInInh=pd(pd(:,4)>0,2); %indices of the inputs in inhibiting links
IdxOutInh=pd(pd(:,4)>0,1); %indices of the outputs in inhibiting links

% % Get Boolean gate indices in for-loop
BoolMax=max(param_index(:,5)); %the number of boolean gates

% store the parameters in 'estim'
estim.Interactions=Interactions;
estim.Input=Input_vector;
estim.Input_idx=Input_index;
estim.Output=Output_vector;
estim.Output_idx=Output_index;

%%% quick fix, search why SD_vector does not have the right length...
if ~isempty(SD_vector)
    SD_vector=SD_vector(:,1:size(Output_vector,2));
else
    SD_vector=zeros(size(Output_vector));
end
estim.SD=SD_vector;
estim.state_names=state_names;
estim.NrStates=size(state_names,2);
estim.NrParams=size(param_vector',2);
estim.param_index=param_index;
estim.param_vector=param_vector;
estim.ma=ma;
estim.mi=mi;
estim.Aeq=Aeq;
estim.A=A;
estim.beq=beq;
estim.b=b;
estim.LB=LB;
estim.UB=UB;
estim.kInd=kInd;
estim.IdxInAct=IdxInAct;
estim.IdxOutAct=IdxOutAct;
estim.IdxInInh=IdxInInh;
estim.IdxOutInh=IdxOutInh;
estim.BoolMax=BoolMax;
estim.BoolIdx=BoolIdx;
estim.BoolOuts=BoolOuts;
estim.FixBool=FixBool;
estim.Annotation = Annotation;
disp('... ... ...')

disp(['Network loaded: ', num2str(estim.NrStates), ' nodes and ', num2str(size(estim.Interactions,1)), ' interactions (', num2str(BoolMax), ' Boolean gates)'])
disp(['Data loaded: ', num2str(size(Input_vector,2)), ' inputs, ', num2str(size(Output_vector,2)), ' outputs and ', num2str(size(Input_vector,1)), ' experimental conditions'])
disp(['The model has ', num2str(estim.NrParams), ' parameters for ', num2str(prod(size(Output_vector))), ' datapoints'])

end
