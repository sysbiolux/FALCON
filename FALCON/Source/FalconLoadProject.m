function FalconLoadProject(Name)

Point = find(ismember(Name,'.'), 1, 'last');
PartName = Name(1:Point-1);
NewName = [PartName, '.mat'];
copyfile(Name, NewName);
evalin('base',['load ', NewName]);
delete(NewName)

end

