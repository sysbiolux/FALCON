function FalconSaveProject(Name)

save(Name);
evalin('base',['save ', Name, '.mat']);
copyfile([Name, '.mat'], [Name, '.falcon'])
delete([Name, '.mat'])

end
