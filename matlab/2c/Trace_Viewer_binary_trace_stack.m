% Convert raw traces to binary traces, and create a stack of traces sorted by
% initial event time

% Modified from Digvijay Singh  ( dgvjay@illinois.edu)
% Edited Ashlee Feng (xfeng17@jhu.edu)
% Last updated July 11, 2021

minInt=100;
ymax=1000;
threshold = 300;

GenericFileType='.traces';   

File_name = 'hel1.traces';
File_id=fopen(File_name,'r');
if File_id == -1
    fprintf(strcat('Error: File ', File_name, ' does not exist.\n'));
    return
end

Timeunit=0.1; % seconds

GammaFactor=1.0;

ChannelLeakage=0.12;

% Extracting important information from .traces binary file
Length_of_the_TimeTraces=fread(File_id,1,'int32');% This fetches the total duration of the imaging that was carried out for the

num_traces=fread(File_id,1,'int16'); 
num_molecules = num_traces / 2;

Raw_Data=fread(File_id,num_traces*Length_of_the_TimeTraces,'int16');
junk_traces_file = fopen('junk.csv');

if junk_traces_file == -1
    junk_traces = {};
else
    junk_traces = textscan(junk_traces_file, '%d', 'Delimiter', ',');
end
disp('Done reading data');
fclose(File_id);

Index_of_SelectedSpots=(1:num_traces*Length_of_the_TimeTraces);
DataMatrix=zeros(num_traces,Length_of_the_TimeTraces);
Donors=zeros(num_traces/2,Length_of_the_TimeTraces);
Acceptors=zeros(num_traces/2,Length_of_the_TimeTraces);
DataMatrix(Index_of_SelectedSpots)=Raw_Data(Index_of_SelectedSpots);

for i=1:(num_traces/2)
    Donors(i,:)=DataMatrix(i*2-1,:);  
    Acceptors(i,:)=GammaFactor.*DataMatrix(i*2,:);
end

BinaryMatrix = Donors < threshold;

TimeSeries=(0:(Length_of_the_TimeTraces-1))*Timeunit;

% Find time when initial event occurs
initial_event = zeros(num_traces/2, 1);

TracesCounter = 0;

while TracesCounter < num_traces/2
    
    TracesCounter = TracesCounter + 1 ;
%     if (find(junk_traces{:} == TracesCounter) > 0)
%         continue
%     end
    bin_trace = BinaryMatrix(TracesCounter, :);
    for t = 1: Length_of_the_TimeTraces
        if bin_trace(t) == 0
            initial_event(TracesCounter) = t;
            break
        end
    end
end

[initial_event, idx] = sort(initial_event, 'descend');
sorted_binary_matrx = BinaryMatrix(idx, :);
n_good = num_traces/2 - size(junk_traces{:}, 1);

figure
imshow(sorted_binary_matrx(1:n_good, :));
