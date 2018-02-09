function [MeasFile] = FalconData2File(estim)
% FalconData2File creates a measurement file from the data in estim.
% The file can then be used in other functions of the Falcon pipeline.

% :: Input values ::
% estim                     complete model definition

% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu


NamesIn=estim.state_names(estim.Input_idx);

NamesOut=estim.state_names(estim.Output_idx);
%Test for presence of excel to select the I/O method
useexcel = isExcelPresent();

delete 'tempMeasFile.xls';
if useexcel    
    xlswrite('tempMeasFile.xls', cellstr(estim.Annotation) , 1, 'A1');
    xlswrite('tempMeasFile.xls', NamesIn , 1, 'B1');
    xlswrite('tempMeasFile.xls', estim.Input, 1, 'B2');
    xlswrite('tempMeasFile.xls', NamesOut, 2, 'A1');
    xlswrite('tempMeasFile.xls', estim.Output, 2, 'A2');
    xlswrite('tempMeasFile.xls', NamesOut, 3, 'A1');
    xlswrite('tempMeasFile.xls', estim.SD, 3, 'A2');
else    
    setupxlwrite();
    xlswrite('tempMeasFile.xls', cellstr(estim.Annotation) , 1, 'A1');
    xlwrite('tempMeasFile.xls', NamesIn , 1, 'B1');
    xlwrite('tempMeasFile.xls', estim.Input, 1, 'B2');
    xlwrite('tempMeasFile.xls', NamesOut, 2, 'A1');
    xlwrite('tempMeasFile.xls', estim.Output, 2, 'A2');
    xlwrite('tempMeasFile.xls', NamesOut, 3, 'A1');
    xlwrite('tempMeasFile.xls', estim.SD, 3, 'A2');
end

MeasFile='tempMeasFile.xls';

end
