function setupxlwrite()
%setupxlwrite setup all folders and pathes necessary to use xlwrite from Matlab file Exchange.
%   setupxlwrite will check whether xlwrite is on the folder, and if not
%   add the FALCONROOT/ThirdParty folder to the path. 
%   it will further check, whether the apache-poi libraries are on the
%   java class path (i.e. if the poi Workbooks classes can be accessed),
%   and if not, it will add the jars in the poi_library folder in
%   FALCONROOT/ThirdParty to the javaclass path.

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

