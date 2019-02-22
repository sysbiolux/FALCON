function [k] = FalconIC(estim, varargin)
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

DefaultMode = 'uniform';
DefaultStd = 0.167;
% to obtain different random starts on fresh Matlab sessions
Clock = clock(); Clock(6) = Clock(6)*1000;
DefaultSeed = sum(Clock);

ExpectedModes = {'uniform', 'normal', 'scratch'};
ValidStd = @(x) isnumeric(x) && isscalar(x) && (x > 0);
% parsing inputs
p = inputParser;
addRequired(p, 'estim');
addOptional(p, 'Mode', DefaultMode,...
        @(x) any(validatestring(x, ExpectedModes)));
addOptional(p, 'Std', DefaultStd, ValidStd);
addOptional(p, 'Seed', DefaultSeed, ValidStd);
parse(p, estim, varargin{:});

Mode = p.Results.Mode;
Std = p.Results.Std;
Seed = p.Results.Seed;

rng(Seed)

if strcmp(Mode, 'uniform') 
    k = rand(1, size(estim.param_vector, 1));
elseif strcmp(Mode, 'normal')
    k = min(max((Std.*randn(1, size(estim.param_vector, 1)) + 0.5), 0), 1);
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