function [k] = FalconIC(varargin)
% FalconIC automatically assigns initial parameter values which are drawn from a uniform or normal distribution. 
% Initial conditions are subsequently scaled according to the equality and inequality constrains.
% k=FalconIC(estim,IC_Dist,Std)
%
% :: Input ::
% estim     complete model definition
% IC_Dist   type of distribution for initial conditions
%           1 = uniform parameter distribution 
%           2 = normal parameter distribtion 
% Std[0-1]  assigned standard deviation for normal distribution
%           (default value = 0.167 which makes the 99%IC=[0-1]
%
% :: Output ::
% k          initial guess for parameter values
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu


estim=varargin{1};

switch nargin
    case 1
        Mode='uniform';
    case 2
        Mode=varargin{2};
        Std=0.167;
    case 3
        Mode=varargin{2};
        Std=varargin{3};
end
    

if strcmp(Mode,'uniform') 

    k=rand(1,size(estim.param_vector,1));
    
elseif strcmp(Mode,'normal')
    
    k=min(max((Std.*randn(1,size(estim.param_vector,1))+0.5),0),1);

end
    
    
% scaling down parameters so that they satisfy the constrains
if ~isempty(estim.A)
    for neq=1:size(estim.A,1)
        IdxConst=find(estim.A(neq,:)>0);
        k(IdxConst)=k(IdxConst)./length(IdxConst);
    end
end
if ~isempty(estim.Aeq)
    for eq=1:size(estim.Aeq,1)
        IdxConst=find(estim.Aeq(eq,:)>0);
        k(IdxConst)=k(IdxConst)./sum(k(IdxConst));
    end
end

end