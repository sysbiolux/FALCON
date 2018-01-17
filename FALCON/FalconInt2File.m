function [] = FalconInt2File(Interactions,fileName)
% still needs a help fuction
%
fileID=fopen(fileName,'w');
formatSpec = '%s\t%s\t%s\t%s\t%s\t%s\n';
[nrows,~] = size(Interactions);
for row = 1:nrows
    fprintf(fileID,formatSpec,Interactions{row,2:end});
end
fclose('all');

end