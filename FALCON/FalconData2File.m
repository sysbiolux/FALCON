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

delete 'tempMeasFile.xlsx';

xlswrite('tempMeasFile.xlsx', NamesIn , 1, 'A1');
xlswrite('tempMeasFile.xlsx', estim.Input, 1, 'A2');
xlswrite('tempMeasFile.xlsx', NamesOut, 2, 'A1');
xlswrite('tempMeasFile.xlsx', estim.Output, 2, 'A2');
xlswrite('tempMeasFile.xlsx', NamesOut, 3, 'A1');
xlswrite('tempMeasFile.xlsx', estim.SD, 3, 'A2');


MeasFile='tempMeasFile.xlsx';

end