function [DidChange, I2]=FalconExpand(Interactions)
% FalconExpand allows for the automatic expansion of logical networks in order to facilitate rapid computation
% [DidChange, I2]=FalconExpand(Interactions)
%
% :: Input values ::
% Interactions      list of interactions, as provided in a InputFile (see FalconMakeModel)
%
% :: Output value(s) ::
% DidChange         toggle switch for model expansion (complete/incomplete)
% I2                list of expanded interactions
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

DidChange=0; %no change... yet
I2=Interactions;
ListOut=unique(Interactions(:,4)); %list of output nodes

for i=1:length(ListOut) %for each output node
    ThisOut=ListOut(i); %name of this output node
    Idx=find(ismember(Interactions(:,4), ThisOut)); %indices where this node appears
    ListParams=Interactions(Idx,5); %collects param names for this output node
    UParams=unique(ListParams); %merge them. If there is only one boolean gate this is a single name
    %%%1) Check if number of incoming interactions is either 1 (can be boolean or not) or only direct interactions
    if length(UParams)>1 % if there are >1 params to this node
        ThisGates=unique(Interactions(Idx,6)); %collect gates for this node
        if ~prod(strcmp(ThisGates,'N')) %if not all are 'N': there is a mix of bool gates and direct, or more than one gate
            DidChange=1; %yes we change
            BoolGates=Idx(find(~ismember(Interactions(Idx,6),'N'))); %indices of boolean gates
            OneGate=BoolGates(1); %take one out
            RelatedParam=Interactions(OneGate,5);
            RelatedCName=cell2mat(Interactions(OneGate,4:6)); %parameter name of that boolean gate
            GateCNames={};
            for ii=1:size(Interactions,1)
                GateCNames=[GateCNames; cell2mat(Interactions(ii,4:6))];
            end
            GateIdx=find(ismember(GateCNames,RelatedCName)); %indices where that gate is found
            GateType=unique(Interactions(GateIdx,6)); % fetch the type of gate
            Type=unique(Interactions(GateIdx,3)); %type of interaction (activation or inhibition)
            NewNode=[char(Interactions(GateIdx(1),2))]; %create new node
            Lim=unique(I2(GateIdx(:),7)); %bound of the parameter
            if strcmp(char(GateType),'O'), Link='_OR_'; elseif strcmp(char(GateType),'A'), Link='_AND_'; end

            for g=2:length(GateIdx)
                ThisIdx=GateIdx(g);
                NewNode=[NewNode,Link,char(Interactions(GateIdx(g),2))];
            end
            I2(GateIdx(:),4)=repmat({NewNode},length(GateIdx),1); %first 2 interactions point to the new node
            I2(GateIdx(:),5)=repmat({['k',NewNode]},length(GateIdx),1); %parameter value is fixed to new param (fixed to 1 later on)
            I2(GateIdx(:),3)=repmat({'->'},length(GateIdx),1); %interaction becomes positive in all cases
            nLine=1+size(I2,1); %number of the line
            I2=[I2;{['i',num2str(nLine)], NewNode, cell2mat(Type), cell2mat(ThisOut), cell2mat(RelatedParam), 'N', cell2mat(Lim)}];
            
            break
        end
    end
    %%%2) check if, in boolean gate, there are more than 1 input
    for p=1:length(UParams) %for each parameter
        
        IdxP=Idx(find(ismember(Interactions(Idx,5), (UParams(p))))); %indices where this param appears
        GateType=unique(Interactions(IdxP,6)); % fetch the type of gate
        if ~strcmp(GateType, 'N') %if boolean gate
            if length(IdxP)>2 %if there are more than 2 inputs
                DidChange=1; %yes we change
                if strcmp(char(GateType),'O'), Link='_OR_'; elseif strcmp(char(GateType),'A'), Link='_AND_'; end
                NewNode=[char(Interactions(IdxP(1),2)),Link,char(Interactions(IdxP(2),2))]; %create new node
                I2(Idx(1:2),4)={NewNode;NewNode}; %first 2 interactions point to the new node
                I2(Idx(1:2),5)={['k',NewNode];['k',NewNode]}; %parameter value is fixed to new param (fixed to 1 later on)
                I2(Idx(1:2),3)={'->';'->'}; %interaction becomes positive anyway.
                nLine=1+size(I2,1); %number of the line
%                 Type=unique(I2(Idx(1:2),3)); %type of interaction (activation or inhibition)
                Lim=unique(I2(Idx(1:2),7)); %bound of the parameter
                %new line
                I2=[I2;{['i',num2str(nLine)], NewNode, '->', cell2mat(ThisOut), cell2mat(UParams(p)), cell2mat(GateType), cell2mat(Lim)}];
                
                break
                
            end
            
        end
        
    end
    
end

if DidChange %check for doublons
    for i=size(I2,1):-1:2
        q=I2(i,[2:6]);
        for ii=i-1:-1:1
           r=I2(ii,[2:6]);
           if sum(strcmp(q,r))==5, I2(i,:)=[]; end
        end
    end
end


end