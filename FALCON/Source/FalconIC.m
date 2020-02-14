function [k] = FalconIC(varargin)
% FalconIC automatically assigns initial parameter values which are drawn from a uniform or normal distribution. 
% Initial conditions are subsequently scaled according to the equality and inequality constrains.
% k=FalconIC(estim,IC_Dist,Std)
%
% :: Input ::
% estim     complete model definition
% Mode(optional)   type of distribution for initial conditions
%           'uniform' = uniform parameter distribution 
%           'normal' = normal parameter distribution
%           'scratch' = near-zero (half-normal) parameter distribution
%
% Std[0-1]  assigned standard deviation for normal distribution
%           (default value = 0.167 which makes the 99%IC=[0-1]
% Seed      Seed for the random number generator, define to get
% reproduceable results.
%
% :: Output ::
% k         initial guess for parameter values
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

% to obtain different random starts on fresh Matlab sessions
Clock = clock(); Clock(6) = Clock(6)*1000;
Seed = sum(Clock);

estim = varargin{1};

if nargin > 1
    Seed = varargin{2};
end

rng(Seed)

Mode = estim.IC_Dist;

if strcmp(Mode, 'uniform') 
    k = rand(1, size(estim.param_vector, 1));
elseif strcmp(Mode, 'normal')
    k = min(max((0.167.*randn(1, size(estim.param_vector, 1)) + 0.5), 0), 1);
elseif strcmp(Mode, 'scratch')
    k = zeros(1, size(estim.param_vector, 1))+(rand(1, size(estim.param_vector, 1))./1000000);
end

% scaling down parameters so that they satisfy the constrains
if ~isempty(estim.A)
    for neq = 1:size(estim.A, 1)
        IdxConst = find(estim.A(neq, :) > 0);
        k(IdxConst) = k(IdxConst)./length(IdxConst);
    end
end
if ~isempty(estim.Aeq)
    for eq = 1:size(estim.Aeq,1)
        IdxConst = find(estim.Aeq(eq,:)>0);
        k(IdxConst) = (k(IdxConst)./sum(k(IdxConst))).*estim.beq(eq);
    end
end

end