function []=FalconShowNetwork(estim, PlotAllBiographs, FinalFolderName)
% FalconShowNetwork uses the 'biograph' environment to display the interaction graphs of the networks.
% Display one graph per condition and display the edges width and the node color intensity according to the weights and node values.
% FalconShowNetwork(estim)
%
% :: Input values ::
% estim             complete model definition
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

IG=(estim.ma'>0)|(estim.mi'>0); % Interaction graph
IDs=estim.state_names; % Names

Net=biograph(IG,IDs); % Make network file

if ~isempty(estim.MeanStateValueAll)
    
    estim.MeanStateValueAll(estim.MeanStateValueAll>1)=1;
    estim.MeanStateValueAll(estim.MeanStateValueAll<0)=0;
    
    EdgeVals=[];
    EdgeTypes=[];
    
    for in=1:length(estim.state_names)
        for out=1:length(estim.state_names)
            if in~=out
                if estim.ma(out,in)>0
                    p1=ismember(estim.Interactions(:,4),estim.state_names(out));
                    p2=ismember(estim.Interactions(:,2),estim.state_names(in));
                    found=find(p1+p2==2);
                    
                    % Only take the edge with the higher weight once both types of interaction are assigned
                    for counter=1:length(found)
                        isAct(counter)=strcmp(estim.Interactions(found(counter),3),'->');
                        K(counter)=estim.Interactions(found(counter),5);
                        if str2num(char(K(counter)))>0
                            Val(counter)=str2num(char(K(counter)));
                        else
                            Val(counter)=estim.bestx(find(ismember(estim.param_vector,K(counter))));
                        end
                        Val_position=find(Val==max(Val));
                        Val=Val(Val_position);
                        isAct=isAct(Val_position);
                    end
                    
                    EdgeVals=[EdgeVals;round(Val*100)/100];
                    EdgeTypes=[EdgeTypes;isAct];
                    
                elseif estim.mi(out,in)>0
                    p1=ismember(estim.Interactions(:,4),estim.state_names(out));
                    p2=ismember(estim.Interactions(:,2),estim.state_names(in));
                    found=find(p1+p2==2);
                    
                    for counter=1:length(found)
                        isAct(counter)=strcmp(estim.Interactions(found(counter),3),'->');
                        K(counter)=estim.Interactions(found(counter),5);
                        if str2num(char(K(counter)))>0
                            Val(counter)=str2num(char(K(counter)));
                        else
                            Val(counter)=estim.bestx(find(ismember(estim.param_vector,K(counter))));
                        end
                        Val_position=find(Val==max(Val));
                        Val=Val(Val_position);
                        isAct=isAct(Val_position);
                    end
                    
                    EdgeVals=[EdgeVals;round(Val*100)/100];
                    EdgeTypes=[EdgeTypes;isAct];
                end
            end
        end
    end
    
    for rep=1:size(estim.MeanStateValueAll,1)
        NodeColors=estim.MeanStateValueAll(rep,:);
        h=view(Net); % Display
        set(h.Nodes,'FontSize',18)
        set(h.Nodes,'Shape','ellipse')
        set(h,'ShowWeights','on')
        set(h,'EdgeFontSize',10)
        for E=1:length(h.Edges)
            set(h.Edges(E),'LineColor',[1-(EdgeTypes(E))/2,(EdgeTypes(E))/1.1,0])
            set(h.Edges(E),'LineWidth',0.5+EdgeVals(E))
            set(h.Edges(E),'Weight',EdgeVals(E))
        end
        for N=1:length(h.Nodes)
            set(h.Nodes(N),'Color',[1,1-NodeColors(N),1])
        end
        
        
        g = biograph.bggui(h);
        f = figure;
        copyobj(g.biograph.hgAxes,f);
        f = get(g.biograph.hgAxes, 'Parent');
        
        if ~PlotAllBiographs
            try                
                print(f, '-dtiff', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.tif'])
                print(f, '-dsvg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.svg'])
                print(f, '-djpeg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.jpg'])

                warning(['Please ignore the error messages above (if any), Biograph figure #' num2str(rep) ' is correctly saved'])
                
                % Close biograph and copied figure after saving
                child_handles = allchild(0);
                names = get(child_handles,'Name');
                BioG = strncmp(['Biograph Viewer'], names, 15);
                close(child_handles(BioG))
                
            catch
                warning('Encountered a graphical memory issue - closing all displayed figures and re-save')
                % Close biograph and copied figure after saving
                child_handles = allchild(0);
                names = get(child_handles,'Name');
                BioG = strncmp(['Biograph Viewer'], names, 15);
                close(child_handles(BioG))
                
                
                try
                    print(f, '-dtiff', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.tif'])
                    print(f, '-dsvg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.svg'])
                    print(f, '-djpeg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.jpg'])

                    warning(['Biograph figure #' num2str(rep) ' is now saved'])
                catch
                    warning(['Not enough graphical memory for saving - Biograph #' num2str(rep) ' is NOT saved'])
                end
            end
            if rep==size(estim.MeanStateValueAll,1)
                rep=1;
                NodeColors=estim.MeanStateValueAll(rep,:);
                h=view(Net); % Display
                set(h.Nodes,'FontSize',18)
                set(h.Nodes,'Shape','ellipse')
                set(h,'ShowWeights','on')
                set(h,'EdgeFontSize',10)
                for E=1:length(h.Edges)
                    set(h.Edges(E),'LineColor',[1-(EdgeTypes(E))/2,(EdgeTypes(E))/1.1,0])
                    set(h.Edges(E),'LineWidth',0.5+EdgeVals(E))
                    set(h.Edges(E),'Weight',EdgeVals(E))
                end
                for N=1:length(h.Nodes)
                    set(h.Nodes(N),'Color',[1,1-NodeColors(N),1])
                end
            end
        else
            try
                print(f, '-dtiff', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.tif'])
                print(f, '-dsvg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.svg'])
                print(f, '-djpeg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.jpg'])
                warning(['Please ignore the error messages above (if any), Biograph figure #' num2str(rep) ' is correctly saved'])
                                
            catch
                warning('Encountered a graphical memory issue - closing all displayed figures and re-save')
                % Close biograph and copied figure after saving
                child_handles = allchild(0);
                names = get(child_handles,'Name');
                BioG = strncmp(['Biograph Viewer'], names, 15);
                close(child_handles(BioG))
                                
                try
                    print(f, '-dtiff', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.tif'])
                    print(f, '-dsvg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.svg'])
                    print(f, '-djpeg', [FinalFolderName, filesep, 'Biograph_Exp' num2str(rep) '.jpg'])
                    warning(['Biograph figure #' num2str(rep) ' is now saved'])
                catch
                    warning(['Not enough graphical memory for saving - Biograph #' num2str(rep) ' is NOT saved'])
                end
            end
        end
    end
    
else
    h=view(Net); % Display
    set(h.Nodes,'FontSize',15)
    
end

end