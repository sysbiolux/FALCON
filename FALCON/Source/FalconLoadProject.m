function FalconLoadProject(Name)

Point = find(ismember(Name,'.'), 1, 'last');
NewName = [Name(1:Point-1), '.mat'];
copyfile(Name, NewName);

load(NewName)
delete(NewName)
end

