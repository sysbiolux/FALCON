function  setupxlwrite()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Check whether xlwrite is on the path. if not, add the path.

if strcmp(which('xlwrite'),'')
    falconfolder = fileparts(which('DriverFalcon'));
    xlwritepath = [falconfolder filesep 'ThirdParty'];
    addpath(xlwritepath);
end
    
    % Check if POI lib is loaded
    
if exist('org.apache.poi.ss.usermodel.WorkbookFactory', 'class') ~= 8 ...
    || exist('org.apache.poi.hssf.usermodel.HSSFWorkbook', 'class') ~= 8 ...
    || exist('org.apache.poi.xssf.usermodel.XSSFWorkbook', 'class') ~= 8
    
    %error('xlWrite:poiLibsNotLoaded',...
    %    'The POI library is not loaded in Matlab.\nCheck that POI jar files are in Matlab Java path!');
    disp('Adding POI Library to Java Path');
    oldfolder = pwd;
    folder = fileparts(which('xlwrite'));
    cd(folder);
    addpath('poi_library');
    javaaddpath(['poi_library' filesep 'poi-3.8-20120326.jar']);
    javaaddpath(['poi_library' filesep 'poi-ooxml-3.8-20120326.jar']);
    javaaddpath(['poi_library' filesep 'poi-ooxml-schemas-3.8-20120326.jar']);
    javaaddpath(['poi_library' filesep 'xmlbeans-2.3.0.jar']);
    javaaddpath(['poi_library' filesep 'dom4j-1.6.1.jar']);
    cd(oldfolder);
end

end

