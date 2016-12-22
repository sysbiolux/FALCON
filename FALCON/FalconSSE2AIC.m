function [AIC] = FalconSSE2AIC(estim, SSE)

N = numel(estim.Output);
np= numel(estim.param_vector);
AIC = N.*log(SSE./N) + 2*np;
AIC=AIC(:)';
end