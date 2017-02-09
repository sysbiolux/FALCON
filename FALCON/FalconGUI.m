function varargout = FalconGUI(varargin)
% Run FALCON from the Graphical User Interface
%
% :: To be set ::
% Load Network      select the network to be optimized
% Load Data         select the experimental data for your network
%                   1) if you select one experimental file, only this file will considered
%                   2) if you select multiple experimental files,
%                      you have to add in a second pop-up window the fixed edges list
% Save in           select a folder to save the optimization results
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

%      FALCONGUI, by itself, creates a new FALCONGUI or raises the existing
%      singleton*.
%
%      H = FALCONGUI returns the handle to a new FALCONGUI or the handle to
%      the existing singleton*.
%
%      FALCONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FALCONGUI.M with the given input arguments.
%
%      FALCONGUI('Property','Value',...) creates a new FALCONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FalconGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FalconGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%


% Last Modified by GUIDE v2.5 07-Feb-2017 09:43:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FalconGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FalconGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% FalconSplash('SplashScreen.jpg',2000)


% --- Executes just before FalconGUI is made visible.
function FalconGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FalconGUI (see VARARGIN)
% Choose default command line output for FalconGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FalconGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FalconGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Loadbutton1.
function Loadbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to Loadbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.InputFile, handles.InputPath]=uigetfile({'*.*';'*.txt';'*.xls';'*.xlsx'},'Select your network file');
set(handles.NetworkFileDisplay,'String',[handles.InputPath, handles.InputFile])
guidata(hObject,handles);

% --- Executes on button press in Loadbutton2.
function Loadbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to Loadbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.MeasFile, handles.DataPath]=uigetfile({'*.*';'*.txt';'*.xls';'*.xlsx'},'Select your data file','MultiSelect','on');
set(handles.DataFileDisplay,'String',[handles.DataPath, handles.MeasFile])
if iscell(handles.MeasFile)
    handles.FixedEdgesList=uigetfile({'*.*';'*.txt';'*.xls';'*.xlsx'},'Select your fixed edges file');
end
if iscell(handles.MeasFile)
    Min_=max(cell2mat((strfind(handles.MeasFile,'_'))'),[],2);
    Max_=max(cell2mat((strfind(handles.MeasFile,'.'))'),[],2);
    handles.ContextsList={};
    for Context=1:length(handles.MeasFile)
        M=char(handles.MeasFile(Context));
        handles.ContextsList=[handles.ContextsList;M(Min_(Context)+1:Max_(Context)-1)];
    end
end
guidata(hObject,handles)

% --- Executes on button press in SaveUI.
function SaveUI_Callback(hObject, eventdata, handles)
% hObject    handle to SaveUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SaveFolderName=uigetdir;
set(handles.SaveAsDisplay,'String',handles.SaveFolderName)
guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function optRoundUI_Callback(hObject, eventdata, handles)
% hObject    handle to optRoundUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function optRoundUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optRoundUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function HLboundUI_Callback(hObject, eventdata, handles)
% hObject    handle to HLboundUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function HLboundUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HLboundUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ParallelUI.
function ParallelUI_Callback(hObject, eventdata, handles)
% hObject    handle to ParallelUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function MaxFunEvalUI_Callback(hObject, eventdata, handles)
% hObject    handle to MaxFunEvalUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function MaxFunEvalUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxFunEvalUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxIterUI_Callback(hObject, eventdata, handles)
% hObject    handle to MaxIterUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function MaxIterUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxIterUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ResultsSingleUI.
function ResultsSingleUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsSingleUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ResultsAllUI.
function ResultsAllUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsAllUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ResultsStatesUI.
function ResultsStatesUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsStatesUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ResultsHeatmapUI.
function ResultsHeatmapUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsHeatmapUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ResultsConvergenceUI.
function ResultsConvergenceUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsConvergenceUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ResultsSummaryUI.
function ResultsSummaryUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsSummaryUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ResultsSummaryUI

% --- Executes on button press in ResamplingUI.
function ResamplingUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResamplingUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in NormalUI.
function NormalUI_Callback(hObject, eventdata, handles)
% hObject    handle to NormalUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LPSAUI.
function LPSAUI_Callback(hObject, eventdata, handles)
% hObject    handle to LPSAUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in DoKOUI.
function DoKOUI_Callback(hObject, eventdata, handles)
% hObject    handle to DoKOUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in DoKONode.
function DoKONode_Callback(hObject, eventdata, handles)
% hObject    handle to DoKONode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ResampleRepUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResampleRepUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function ResampleRepUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ResampleRepUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in AllPlots.
function AllPlots_Callback(hObject, eventdata, handles)
% hObject    handle to AllPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AllPlots

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LPSARepUI_Callback(hObject, eventdata, handles)
% hObject    handle to LPSARepUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function LPSARepUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LPSARepUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ResultsfminconConvergenceUI.
function ResultsfminconConvergenceUI_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsfminconConvergenceUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LPSA_fast.
function LPSA_fast_Callback(hObject, eventdata, handles)
% hObject    handle to LPSA_fast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LPSA_fast


% --- Executes on button press in BiographUI.
function BiographUI_Callback(hObject, eventdata, handles)
% hObject    handle to BiographUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BiographUI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in AnalyseUI.
function AnalyseUI_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyseUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.AbortUI,'Value',0)

% FalconInstall; clc % In case the toolbox is not yet installed

optRound=str2double(get(handles.optRoundUI,'String'));
HLbound=str2double(get(handles.HLboundUI,'String'));
% SSthresh=str2double(get(handles.SSthreshUI,'String'));
SSthresh=eps;
Parallel=get(handles.ParallelUI,'Value');
MaxFunEval=str2double(get(handles.MaxFunEvalUI,'String'));
MaxIter=str2double(get(handles.MaxIterUI,'String'));
ResultsForm(1)=get(handles.ResultsSingleUI,'Value');
ResultsForm(2)=get(handles.ResultsAllUI,'Value');
ResultsForm(4)=get(handles.ResultsStatesUI,'Value');
ResultsForm(3)=get(handles.ResultsHeatmapUI,'Value');
ResultsForm(5)=get(handles.ResultsConvergenceUI,'Value');
ResultsSummary=get(handles.ResultsSummaryUI,'Value');
FminconConv=get(handles.ResultsfminconConvergenceUI,'Value');
DoResample=get(handles.ResamplingUI,'Value');
ResampleRep=str2double(get(handles.ResampleRepUI,'String'));
UseNormal=get(handles.NormalUI,'Value');
DoLPSA=get(handles.LPSAUI,'Value');
LPSARep=str2double(get(handles.LPSARepUI,'String'));
DoKO=get(handles.DoKOUI,'Value');
DoKONo=get(handles.DoKONode,'Value');
UseFast=get(handles.LPSA_fast,'Value');
DoBiograph=get(handles.BiographUI,'Value');
AllPlots=get(handles.AllPlots,'Value');

%%%%%%%

%%% set-up
set(handles.ProgressDisplay2,'String','Creating Model...'); drawnow
set(handles.ProgressDisplay,'String',''); drawnow
if ischar(handles.MeasFile)    
    estim=FalconMakeModel([handles.InputPath filesep handles.InputFile],[handles.DataPath filesep handles.MeasFile],HLbound,1); %make the model
elseif iscell(handles.MeasFile)
    MeasFileList=handles.MeasFile;
    [estim, stamp] = FalconMakeGlobalModel(handles.InputFile,handles.FixedEdgesList,MeasFileList,handles.ContextsList,HLbound,1);
    MeasFile=['Results_' stamp '_.xlsx'];
end

estim.options = optimoptions('fmincon','TolCon',1e-6,'TolFun',1e-6,'TolX',1e-10,'MaxFunEvals',MaxFunEval,'MaxIter',MaxIter); % Default
estim.SSthresh=SSthresh;
if UseNormal
    UN='normal';
else
    UN='uniform';
end

%%% optimization
toc_all=[];
x_all=[];
fval_all=[];

StopCommand=0;

if ~StopCommand
    set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
    set(handles.ProgressDisplay2,'String','Optimizing...'); drawnow
    if Parallel
        parfor counter=1:optRound
            set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
            set(handles.ProgressDisplay2,'String',['Optimization ',num2str(counter),' out of ',num2str(optRound)]); drawnow
            tic
            k=FalconIC(estim,UN);
            [xval,fval]=FalconObjFun(estim,k);
            Fitting_Time=toc;
            toc_all=[toc_all; Fitting_Time];
            x_all=[x_all; xval];
            fval_all=[fval_all; fval];
            
        end
    else
        for counter=1:optRound
            set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
            set(handles.ProgressDisplay2,'String',['Optimization ',num2str(counter),' out of ',num2str(optRound)]); drawnow
            tic
            k=FalconIC(estim,UN);
            [xval,fval]=FalconObjFun(estim,k);
            Fitting_Time=toc;
            toc_all=[toc_all; Fitting_Time];
            x_all=[x_all; xval];
            fval_all=[fval_all; fval];
            
        end

    end
    fxt_all=[fval_all x_all toc_all];
    beep; pause(0.5); beep;
else
    error('Terminated by the user')
end

estim.MaxTime = mean(fxt_all(:,end))*3;
estim.Results.Optimisation.FittingCost = fxt_all(:,1);
estim.Results.Optimisation.FittingTime = fxt_all(:,end);
estim.Results.Optimisation.ParamNames = estim.param_vector;
estim.Results.Optimisation.StateNames = estim.state_names;

%%% Retrieving the results
StopCommand=0;

if ~StopCommand
    set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
    set(handles.ProgressDisplay2,'String',['Optimization over. Computing results...']); drawnow
    [bestx,meanx,stdx]=FalconResults(fxt_all,estim.param_vector);
    estim.Results.Optimisation.BestParams = bestx;
    
    warning off
    if ResultsSummary
        ResultFileName=['Summary',datestr(now, 'yyyymmddTHHMMSS'),'.xlsx'];
        useexcel = isExcelPresent();
        if useexcel
            %Try the Excel Backend,             
            xlswrite([handles.SaveFolderName, filesep, ResultFileName],{'Parameter','Best','Average','Std'},'1','A1');
            xlswrite([handles.SaveFolderName, filesep, ResultFileName],estim.param_vector,'1','A2');
            xlswrite([handles.SaveFolderName, filesep, ResultFileName],[bestx',meanx',stdx'],'1','B2');    
        else
            %Otherwise save as a csv
            ResultFileName=['Summary',datestr(now, 'yyyymmddTHHMMSS'),'.csv'];            
            tab = table(estim.param_vector,bestx',meanx',stdx','VariableNames',{'Parameter','Best','Average','Std'});
            writetable(tab,[handles.SaveFolderName, filesep, ResultFileName],'Delimiter',',');
        end
    end
    warning on
else
    error('Terminated by the user')
end

%%% Fitting cost evolution

StopCommand=0;

if ~StopCommand
    if FminconConv
        FalconFitEvol(estim,UN,handles.SaveFolderName);
    end
else
    error('Terminated by the user')
end

%%% Analyzing the results

StopCommand=0;

if ~StopCommand
    set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
    set(handles.ProgressDisplay2,'String',['Compiling results...']); drawnow
    [MeanStateValueAll, StdStateValueAll, MeanCostAll, StdCostAll, estim] = FalconSimul(estim,bestx,ResultsForm,handles.SaveFolderName);
    estim.MeanStateValueAll=MeanStateValueAll; estim.bestx=bestx;
else
    error('Terminated by the user')
end

%%% Biograph

StopCommand=0;

if ~StopCommand
    if DoBiograph
        set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
        set(handles.ProgressDisplay2,'String',['Plotting Network structures with weighted edges...']); drawnow
        FalconShowNetwork(estim, AllPlots, handles.SaveFolderName)
    end
else
    error('Terminated by the user')
end

%%% Obtaining robust estimates for the parameters by resampling

StopCommand=0;

if ~StopCommand
    if DoResample
        set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
        set(handles.ProgressDisplay2,'String',['Computing Resampling of datapoints...']); drawnow
        handles.SaveFolderName
        CV_cutoff = 10; % Percent cut-off for large coefficient of variation
        [~,~,estim]=FalconResample(estim, bestx, 3, ResampleRep, CV_cutoff, handles.SaveFolderName);
                
    end
else
    error('Terminated by the user')
end

%%% Identifiability Analysis (LPSA)

StopCommand=0;

if UseFast
    UseFast='fast';
else
    UseFast='complete';
end

if ~StopCommand
    if DoLPSA
        set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
        set(handles.ProgressDisplay2,'String',['Computing identifiability of parameters...']); drawnow
        [~,estim]=FalconLPSA(estim, bestx, handles.MeasFile, HLbound,1,LPSARep,UseFast,Parallel,handles.SaveFolderName);
        
    end
else
    error('Terminated by the user')
end

%%% Knock-Outs

StopCommand=0;

if ~StopCommand
    if DoKO
        set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
        set(handles.ProgressDisplay2,'String',['Computing Systematic Knock-Outs...']); drawnow
        estim=FalconKO(estim, bestx, fxt_all,handles.MeasFile,HLbound,1,handles.SaveFolderName);
    end
else
    error('Terminated by the user')
end

StopCommand=0;

if ~StopCommand
    if DoKONo
        set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
        set(handles.ProgressDisplay2,'String',['Computing Systematic Node Knock-Outs...']); drawnow
        estim=FalconKONodes(estim, bestx, fxt_all,handles.MeasFile,HLbound,1,handles.SaveFolderName);
            end
else
    error('Terminated by the user')
end


%%% Stop button interpretation 
if ~StopCommand
    set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
    set(handles.ProgressDisplay2,'String',['Computations over and everything is all right!']); drawnow
else
    set(handles.ProgressDisplay,'String',get(handles.ProgressDisplay2,'String')); drawnow
    set(handles.ProgressDisplay2,'String',['Computations aborted by user!']); drawnow
    error('Terminated by the user')
end

disp(' ')
disp('================================================')
disp('The optimisation and analyses are completed')
disp('================================================')

disp(' ')
disp('Please also check the results in "estim.Results"')
disp(estim.Results)
warning off
save([handles.SaveFolderName filesep 'estim_Results'])
warning on
disp('================================================')
