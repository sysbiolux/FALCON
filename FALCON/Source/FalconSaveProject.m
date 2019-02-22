function FalconSaveProject(Name)

save(Name);
copyfile([Name, '.mat'], [Name, '.falcon'])
delete([Name])

end
