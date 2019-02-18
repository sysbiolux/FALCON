function [estim2]=FalconPreProcess(estim, varargin)
% FalconPreProcess transforms the data in the estim structure. Use it to
% normalize min-max, bootstrap in any or both dimensions, randomize any
% fraction of the data, and sub- or over-sample.
% FalconPreProcess(estim,[Option, Value])
% 
% :: Input values ::
% estim             complete model definition
% Option           : type of distribution for the initial condition
% Value            : argumrnt for the option
%
% :: Output values ::
% estim             updated model definition
% Examples:
% [estim2]=FalconPreProcess(estim, 'normalize', [0 1])
% [estim2]=FalconPreProcess(estim, 'bootstrap', 'rows')
% [estim2]=FalconPreProcess(estim, 'bootstrap', 'columns')
% [estim2]=FalconPreProcess(estim, 'bootstrap', 'both')
% [estim2]=FalconPreProcess(estim, 'randomize', 0.5)
% [estim2]=FalconPreProcess(estim, 'subsample', 0.5)
% [estim2]=FalconPreProcess(estim, 'subsample', 2.5)
% [estim2]=FalconPreProcess(estim, 'normalize', [-0.2 0.9], 'bootstrap', 'both','randomize', 0.5, 'subsample', 0.5)
% [estim2]=FalconPreProcess(estim, 'noise', 1)
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu
% FALCON July 2017

Norm=0;
BS_row=0;
BS_col=0;
Rand=0;
Subs=0;
Noise=0;

estim2=estim;

for c=1:numel(varargin)
    Option=varargin{c};
    
    if strcmp('normalize',Option)
        Norm=1; Range=varargin{c+1};
        t_min=min(Range); t_max=max(Range);

    elseif strcmp('bootstrap',Option)
        if strcmp(varargin{c+1},'rows')
            BS_row=1;
        elseif strcmp(varargin{c+1}, 'columns')
            BS_col=1;
        elseif strcmp(varargin{c+1}, 'both')
            BS_row=1; BS_col=1;
        end

    elseif strcmp('randomize', Option)
        Rand=1;
        FractR=varargin{c+1};

    elseif strcmp('subsample', Option)
        Subs=1;
        FractS=varargin{c+1};
    elseif strcmp('noise', Option)
        Noise=1;
        FractN=varargin{c+1};
    end

end

X=estim.Input;
Y=estim.Output;

if Norm
    disp('normalizing...')
    Idx=(std(X,0,1))~=0;
    X(:,Idx)=(X(:,Idx)-min(X(:,Idx)))./(max(X(:,Idx))-min(X(:,Idx)));
    X=(X.*(t_max-t_min))+t_min;
    Y=(Y-min(Y))./(max(Y)-min(Y));
    Y=(Y.*(t_max-t_min))+t_min;
    estim2.Input=X;
    estim2.Output=Y;
    %%% TODO: still need to normalize SD
end

[nRows,nColsX]=size(X);
[~     ,nColsY]=size(Y);

if BS_row
    disp('bootstrapping rows...')
    Chosen=ceil(rand(1,nRows).*nRows);
    estim2.Output=Y(Chosen,:);
    estim2.SD=estim.SD(Chosen,:);
end

if BS_col
    disp('bootstrapping columns...')
    Chosen=ceil(rand(1,nColsY).*nColsY);
    estim2.Output=Y(:,Chosen);
    estim2.SD=estim.SD(Chosen,:);
end

if Rand
    disp('randomizing...')
    Perm_Rows=randperm(nRows);
    Perm_ColsY=randperm(nColsY);
    
    Y2=Y(Perm_Rows, Perm_ColsY);
    SD2=estim.SD(Perm_Rows, Perm_ColsY);
    
    My=rand(nRows,nColsY)<FractR;
    
    Y(My)=Y2(My);
    estim2.SD(My)=SD2(My);
    estim2.Output=Y;
end

if Subs
    disp('subsampling...')
    Chosen=randperm(nRows);
    Chosen=repmat(Chosen,1,ceil(FractS));
    Chosen=Chosen(1:ceil(FractS*nRows));
    X=X(Chosen,:);
    Y=Y(Chosen,:);
    estim2.Input=X;
    estim2.Output=Y;
    estim2.Input_idx=estim.Input_idx(Chosen,:);
    estim2.Output_idx=estim.Output_idx(Chosen,:);
    estim2.SD=estim.SD(Chosen,:);
    estim2.Annotation=estim.Annotation(Chosen,:);
end

if Noise
    disp('...adding noise...')
    Y=estim.Output;
    SD=estim2.SD;
    Y2=randn(size(Y)).*SD.*FractN;
    Y=max(min((Y+Y2),1),0);
    estim2.Output=Y;
end



end