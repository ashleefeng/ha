Directory_of_Traces = input('Directory:', 's');
if isempty(Directory_of_Traces)
    Directory_of_Traces = pwd;
end

cd(Directory_of_Traces);

list_of_files = dir([fullfile(Directory_of_Traces) '\\*.pks']);
len = length(list_of_files);
num_pks = zeros(len, 1);

for file_num = 1:len
    fileID = fopen(list_of_files(file_num).name, 'r');
    
    for i = 1:6
        temp_read = fgets(fileID);
    end
    
    num_pks(file_num) = str2double(temp_read);
    fclose(fileID);
end

avg = mean(num_pks);
stdev = std(num_pks);

fprintf('Spot count: %f +- %f\n', avg, stdev);


