function useexcel = isExcelPresent()
%isExcelPresent checks, whether there is a working excel instance present on this machine
%   The function uses the matlab function to obtain an excel instance, and,
%   if that fails returns false, otherwise true
% isExcelPresent() - return true if an excel instance can be obtained,
% false otherwise

useexcel = true;
try
    Excel = matlab.io.internal.getExcelInstance; %This fails if no excel instance exists.        
catch exc   
    useexcel = false;
end

end

