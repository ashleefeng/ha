% Written by Digvijay Singh  ( dgvjay@illinois.edu)
% Edited by Ashlee Feng (xfeng17@jhu.edu)

clear
warning('off','all')

fclose('all');
minInt=100;
ymax=800;

Directory_of_TracesFiles=input('Directory: ','s');

if isempty(Directory_of_TracesFiles)
    Directory_of_TracesFiles=pwd;
end

cd(Directory_of_TracesFiles);

GenericFileType='.traces';

FileIndexNumber=input('Index: ');

if isempty(FileIndexNumber)
    FileIndexNumber=1;
end

FileIndexNumber=num2str(FileIndexNumber);

TheFileThatWillbeAnalysed = strcat('hel', FileIndexNumber, GenericFileType);
fprintf(strcat(TheFileThatWillbeAnalysed, '\n'));

File_name = strcat('hel', FileIndexNumber, GenericFileType);
File_id=fopen(File_name,'r');
if File_id == -1
    fprintf(strcat('Error: File ', File_name, ' does not exist.\n'));
    return
end

Timeunit=input('Enter the value of the time unit i.e. frame rate [Default=0.1 sec] ');

if isempty(Timeunit)
    Timeunit=0.1;
end

GammaFactor=input('Enter the value of the Gamma Factor [Default is 1.0] :  ');

if isempty(GammaFactor)
    GammaFactor=1.0;
end

ChannelLeakage=input('Enter the value of the Channel Leakage [Default is 0.12] ');
if isempty(ChannelLeakage)
    ChannelLeakage=0.12;
    % T70S: 0.175
end


Length_of_the_TimeTraces=fread(File_id,1,'int32');% This fetches the total duration of the imaging that was carried out for the

disp('The length of the time traces is: ')
disp(Length_of_the_TimeTraces);

num_traces=fread(File_id,1,'int16');

disp('The number of traces in this file is:')
num_molecules = num_traces / 2;
disp(num_molecules);

Raw_Data=fread(File_id,num_traces*Length_of_the_TimeTraces,'int16');
disp('Done reading data');
fclose(File_id);

Index_of_SelectedSpots=(1:num_traces*Length_of_the_TimeTraces);
DataMatrix=zeros(num_traces,Length_of_the_TimeTraces);
Donors=zeros(num_traces/2,Length_of_the_TimeTraces);
Acceptors=zeros(num_traces/2,Length_of_the_TimeTraces);
DataMatrix(Index_of_SelectedSpots)=Raw_Data(Index_of_SelectedSpots);

for i=1:(num_traces/2)
    Donors(i,:)=DataMatrix(i*2-1,:);   %So this will be a matrix where each column will be the Donor time series of each selected spot of the movie
    Acceptors(i,:)=GammaFactor.*DataMatrix(i*2,:); %So this will be a matrix where each column will be the Acceptor time series of each selected spot of the movie
end

TimeSeries=(0:(Length_of_the_TimeTraces-1))*Timeunit;
smoothed_fret_x = zeros(floor(Length_of_the_TimeTraces/3));
for i = 1: floor(Length_of_the_TimeTraces/3)
    smoothed_fret_x(i) = TimeSeries(i*3);
end

pks_fname = strcat('hel', FileIndexNumber, '.pks');
pks_id=fopen(pks_fname,'r');

if pks_id ~= -1
    for i = 1:6
        tline = fgets(pks_id);
    end
    A = fscanf(pks_id,'%f %f %f %f %f, %f %f %f %f %f',[10 Inf]);
    A = A';
    fclose(pks_id);
else
    fprintf(strcat('\nError: ', pks_fname, ' does not exist.\n\n'));
end


newfolder = [FileIndexNumber ' selected traces'];
mkdir(newfolder);

TracesCounter=77;
figure;
hd13=gcf;

counter = 0;


step_count_list = zeros(num_molecules, 1);


tif_fname = strcat('hel', FileIndexNumber, '_ave.tif');
tif_fileID = fopen(tif_fname);

if tif_fileID ~= -1
    fclose(tif_fileID);
else
    fprintf(strcat('\nError: ', tif_fname, ' does not exist.\n\n'));
end

while TracesCounter < num_traces/2
    
    %close all;
    TracesCounter = TracesCounter + 1 ;
    figure(hd13);
    %    ax1=subplot(3,3,[1 2 3]);
%     ax1 = subplot(2, 1, 1);
    ax1=subplot(2,4,[1 2 3]);
    plot(TimeSeries,Donors(TracesCounter,:),'g',...
        TimeSeries,(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:)),'r');
    
    temp=axis;
    temp(2)=TimeSeries(end);
    temp(3)=-100;
    %    temp(4)=ymax; %adjust max y-axis
    axis(temp);
    
    
    xlabel('Time(s)');
    ylabel('Intensity (a.u.)');
    xlim([0 Length_of_the_TimeTraces * Timeunit]);
    TitleNameForThePlot=sprintf('Molecule %d / %d of %s',TracesCounter,num_traces/2,TheFileThatWillbeAnalysed);
    title(TitleNameForThePlot); 
    
%     ax2 = subplot(2, 1, 2);
    ax2=subplot(2, 4,[5 6 7]);
    FRET_Time_Series=(Acceptors(TracesCounter,:)...
        -ChannelLeakage*Donors(TracesCounter,:))...
        ./(Acceptors(TracesCounter,:)...
        -ChannelLeakage*Donors(TracesCounter,:)...
        +(Donors(TracesCounter,:)));
    
    plot(TimeSeries,FRET_Time_Series,'LineWidth',1.2,'Color','b');
    xlabel('Time(s)');
    ylabel('FRET Efficiency');
    xlim([0 Length_of_the_TimeTraces * Timeunit]);
    ylim([-0.2 1.2]);
    temp=axis;temp(2)=TimeSeries(end);
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    
    linkaxes([ax1,ax2],'x');
    
    green = uint8([0 255 0]);

    if tif_fileID ~= -1
        
        [X,map] = imread(['hel' FileIndexNumber '_ave.tif']);
        
        if ~isempty(map)
            Im = ind2rgb(X,map);
        end
        
        dspotx = A(TracesCounter,1); %xcoord
        dspoty = A(TracesCounter,3); %ycoord
        aspotx = A(TracesCounter,6); %xcoord
        aspoty = A(TracesCounter,8); %ycoord
        circleradius=4.5; %dimension
        
        circles = int32([dspotx dspoty circleradius;aspotx aspoty circleradius]);
        
        donorInserter = insertShape(Im, 'circle', [dspotx dspoty circleradius], 'Color', 'cyan');
        acceptorInserter = insertShape(Im, 'circle', [aspotx aspoty circleradius], 'Color', 'cyan');
        
        sz=512;
        subplot(2, 4, 8);
        imshow(imresize(acceptorInserter(int16(max(1,aspoty-20)):int16(min(sz,aspoty+20)),int16(max(1,aspotx-20)):int16(min(512,aspotx+20)),:),4,'nearest'));
        title('Acceptor');
        subplot(2, 4, 4);
        imshow(imresize(donorInserter(int16(max(1,dspoty-20)):int16(min(sz,dspoty+20)),int16(max(1,dspotx-20)):int16(min(512,dspotx+20)),:),4,'nearest'));
        title('Donor');
        zoom on;
        
    end
    
    step_count = input('Category for this trace: 1-HFsteady, 2-HFdynamic, 3-MFdynamic, 4-MFsteady, 5-other, 9-smooth, 0-back\n');
    
    if step_count == 0
        TracesCounter = TracesCounter - 2;
        continue
    end
    
    if step_count == 9
        
        smoothed_fret_y = zeros(floor(Length_of_the_TimeTraces/3));
        
        for i = 1: floor(Length_of_the_TimeTraces/3)
            smoothed_fret_y(i) = mean(FRET_Time_Series(i*3-2:i*3));
        end
        
        ax2=subplot(2, 4,[5 6 7]);
        plot(TimeSeries,FRET_Time_Series,'LineWidth',1.2,'Color',[153/255, 153/255, 153/255]);
        hold on
        plot(smoothed_fret_x,smoothed_fret_y,'LineWidth',1.5,'Color','b');
        hold off
        
        xlabel('Time(s)');
        ylabel('FRET Efficiency');
        xlim([0 Length_of_the_TimeTraces * Timeunit]);
        ylim([-0.1 1.1]);
        
        done_smooth = input('Enter to move on to next trace, 9 to undo\n');
        if done_smooth == 9
            TracesCounter = TracesCounter - 1;
        end
        continue
    end
    
    if ~isempty(step_count)
        step_count_list(TracesCounter) = step_count;
    end
    
end

instr = input('Do you want to analyze trace categories? (Y/N)\n', 's');

if isempty(instr)
    instr = 'Y';
end

if strcmp(instr, 'Y') == 1
    
    figure;
    % create a histogram for photobleaching steps
    count_hist = histogram(step_count_list);
    ch_values = count_hist.Values();
    max_count = size(ch_values, 2);
    xlabel('Category');
    ylabel('Number of Traces');
    
    fname = strcat('hel', FileIndexNumber, '_hfds');
%     png_fname = strcat(fname,'.png');
%     saveas(count_hist, png_fname);
    

    % print the histogram results
    for i_count = 1: max_count
        fprintf('number of category %d traces: %d\n', ...
            i_count-1, ch_values(i_count));
    end
    
    % save the counts as a .csv file
    indices = 1:1:num_molecules;
    print_count_list = [indices' step_count_list];
    csv_fname = strcat(fname, '.csv');
    csvwrite(csv_fname, print_count_list);
    
end

fprintf('Done.\n');

close all;
