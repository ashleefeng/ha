% Written by Digvijay Singh  ( dgvjay@illinois.edu)
% Edited by Olivia Yang (oyang1@jhu.edu) and Ashlee Feng (xfeng17@jhu.edu)
% Last updated Oct 10, 2019

close all
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
TheFileThatWillbeAnalysed = strcat('hel', FileIndexNumber, GenericFileType); %displays the .traces file that will be used to show
fprintf(strcat(TheFileThatWillbeAnalysed, '\n'));

%FilePointer
File_name = strcat('hel', FileIndexNumber, GenericFileType);
File_id=fopen(File_name,'r');
if File_id == -1
    fprintf(strcat('Error: File ', File_name, ' does not exist.\n'));
    return
end

% Define time unit (seconds)
Timeunit=input('Enter the value of the time unit i.e. frame rate [Default=0.1 sec] ');
if isempty(Timeunit)
    Timeunit=0.1;
end

GammaFactor=input('Enter the value of the Gamma Factor [Default is 1.0] :  ');
if isempty(GammaFactor)
    GammaFactor=1.0;
end

% Channel Leakage Value
ChannelLeakage=input('Enter the value of the Channel Leakage [Default is 0.12] ');
if isempty(ChannelLeakage)
    ChannelLeakage=0.12;
    % T70S: 0.175
end

% Extracting important information from .traces binary file
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
    Donors(i,:)=DataMatrix(i*2-1,:);  
    Acceptors(i,:)=GammaFactor.*DataMatrix(i*2,:);
end

TimeSeries=(0:(Length_of_the_TimeTraces-1))*Timeunit;

smoothed_fret_x = zeros(floor(Length_of_the_TimeTraces/3));
for i = 1: floor(Length_of_the_TimeTraces/3)
    smoothed_fret_x(i) = TimeSeries(i*3);
end

pks_fname = strcat('hel', FileIndexNumber, '.pks');
pks_id=fopen(pks_fname,'r');

if pks_id ~= -1
%     for i = 1:6
%         tline = fgets(pks_id);
%     end
    % A = fscanf(pks_id,'%f %f %f %f %f, %f %f %f %f %f',[10 Inf]);
    A = fscanf(pks_id,'%f %f %f %f',[4 Inf]);
    A = A';
    fclose(pks_id);
else
    fprintf(strcat('\nError: ', pks_fname, ' does not exist.\n\n'));
end

newfolder = [FileIndexNumber ' selected traces'];
mkdir(newfolder);

TracesCounter=0;
figure;
hd13=gcf;

DT1=[];DT2=[];DT3=[];
DT1a=[];DT2a=[];DT3a=[];
DT1d=[];DT2d=[];DT3d=[];
DT1f=[];DT2f=[];DT3f=[];

counter = 0;
start_pt = 1;
end_pt = Length_of_the_TimeTraces;

step_count_list = zeros(num_molecules, 1);

tif_fname = strcat('hel', FileIndexNumber, '_ave.tif');
tif_fileID = fopen(tif_fname);

if tif_fileID ~= -1
    fclose(tif_fileID);
else
    fprintf(strcat('\nError: ', tif_fname, ' does not exist.\n\n'));
end

while TracesCounter < num_traces/2
    TracesCounter = TracesCounter + 1 ;
    figure(hd13);
    ax1=subplot(3,3,[1 2 3]);
    plot(TimeSeries,Donors(TracesCounter,:),'g',...
        TimeSeries,(Acceptors(TracesCounter,:) - ChannelLeakage*Donors(TracesCounter,:)),'r');
    hold on
    total = plot(TimeSeries,Acceptors(TracesCounter,:) + (1 - ChannelLeakage) * Donors(TracesCounter,:) + 200, 'k');
    hold off
    
    temp=axis;
    temp(2)=TimeSeries(end);
    temp(3)=-100;
    temp(4)=ymax; %adjust max y-axis
    axis(temp);
    
    
    xlabel('Time(s)');
    ylabel('Intensity (a.u.)');
    xlim([0 Length_of_the_TimeTraces * Timeunit]);
    TitleNameForThePlot=sprintf('Molecule %d / %d of %s',TracesCounter,num_traces/2,TheFileThatWillbeAnalysed);
    title(TitleNameForThePlot);
    
    ax2=subplot(3, 3,[4 5 6]);
    FRET_Time_Series=(Acceptors(TracesCounter,:)...
        -ChannelLeakage*Donors(TracesCounter,:))...
        ./(Acceptors(TracesCounter,:)...
        -ChannelLeakage*Donors(TracesCounter,:)...
        +(Donors(TracesCounter,:)));
    
    plot(TimeSeries(start_pt:end_pt),FRET_Time_Series(start_pt:end_pt),'LineWidth',1.2,'Color','b');
    xlabel('Time(s)');
    ylabel('FRET Efficiency');
    xlim([0 Length_of_the_TimeTraces * Timeunit]);
    ylim([-0.2 1.2]);
    temp=axis;temp(2)=TimeSeries(end);
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    
    linkaxes([ax1,ax2],'x');
    
    subplot(3,3,7);
    histogram(FRET_Time_Series(start_pt:end_pt), -0.2:0.05:1.2);
    xlim([-0.2 1.2]);
    xlabel('FRET Efficiency'); 
    ylabel('Count'); 
    title('Trace FRET histogram');

    % start shape drawing
    
    green = uint8([0 255 0]);
    
    % plot original spot image
    
    if tif_fileID ~= -1
        
        [X,map] = imread(['hel' FileIndexNumber '_ave.tif']);
        
        if ~isempty(map)
            Im = ind2rgb(X,map);
        end
        
%         
%         dspotx = A(TracesCounter,1); %xcoord
%         dspoty = A(TracesCounter,3); %ycoord
%         aspotx = A(TracesCounter,6); %xcoord
%         aspoty = A(TracesCounter,8); %ycoord
        dspotx = A(TracesCounter*2-1,2); %xcoord
        dspoty = A(TracesCounter*2-1,3); %ycoord
        aspotx = A(TracesCounter*2,2); %xcoord
        aspoty = A(TracesCounter*2,3); %ycoord
        circleradius=4.5; %dimension
        
        circles = int32([dspotx dspoty circleradius;aspotx aspoty circleradius]);
        
        donorInserter = insertShape(Im, 'circle', [dspotx dspoty circleradius], 'Color', 'cyan');
        acceptorInserter = insertShape(Im, 'circle', [aspotx aspoty circleradius], 'Color', 'cyan');
        
        sz=512;
        subplot(3, 3, 9);
        imshow(imresize(acceptorInserter(int16(max(1,aspoty-20)):int16(min(sz,aspoty+20)),int16(max(1,aspotx-20)):int16(min(512,aspotx+20)),:),4,'nearest'));
        title('Acceptor');
        subplot(3, 3, 8);
        imshow(imresize(donorInserter(int16(max(1,dspoty-20)):int16(min(sz,dspoty+20)),int16(max(1,dspotx-20)):int16(min(512,dspotx+20)),:),4,'nearest'));
        title('Donor');
        zoom on;
        
    end
    
    choice = input('Press 0 to go back one trace\nPress 1 to save current trace\nPress 2 to subtract background (2 clicks)\nPress 3 to cut out\nPress 4 to save the trace as a picture\nPress 5 to skip to data analysis and exit\nPress 6 to collect dwell times\nPress 7 to choose a trace\nPress 8 to select range for single trace histogram\nPress 9 to smoothen the trace\nPress Enter to move on\n');
    
    if choice==0
        TracesCounter=TracesCounter - 2;
        continue;
    end

    if choice==1
        % use '\' for PC
        fname1=[Directory_of_TracesFiles '/' newfolder  '/hel' FileIndexNumber ' tr' num2str(TracesCounter) '.dat'];
        Output_to_be_saved_inaFile=[TimeSeries()' Donors(TracesCounter,:)' (Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:))'];
        save(fname1,'Output_to_be_saved_inaFile','-ascii') ;
    end
    
    if choice==2
        [x,y]=ginput(2);
        st = floor((x(1)-TimeSeries(1))/Timeunit);
        en = floor((x(2)-TimeSeries(1))/Timeunit);
        i3bg=mean(Donors(TracesCounter,st:en));
        i5bg=mean(Acceptors(TracesCounter,st:en));
        
        Donors(TracesCounter,:) = Donors(TracesCounter,:)-i3bg;
        Acceptors(TracesCounter,:) = Acceptors(TracesCounter,:)-i5bg;
        
        TracesCounter = TracesCounter - 1;
    end
    
    if choice == 3
        TheFileThatWillbeAnalysed_truncated=TheFileThatWillbeAnalysed(1:end-7);
        fprintf('Click twice to select the range to extract\n');
        Done_Cutting_Choice = 0;
        Cut_Counter = 0;
        while Done_Cutting_Choice == 0
            Cut_Counter = Cut_Counter+1;
            [x, y, button] = ginput(2); 
            x(1) = round(x(1)/Timeunit); 
            x(2) = round(x(2)/Timeunit);
            if button(1) ~= button(2)
                fprintf('You should use the same buttons to select each range. Try again.\n')
                continue;
            end
            if button(1) == 1
                % left mouse button
                trace_type = 'type1';
            elseif button(1) == 2
                % middle mouse button
                trace_type = 'type2';
            elseif button(1) == 3
                % right mouse button
                trace_type = 'type3';
            end
            
            %fname1=[Directory_of_TracesFiles '/' newfolder  '/' trace_type '/hel' FileIndexNumber ' tr' num2str(TracesCounter) '-' num2str(Cut_Counter) '.dat'];
            fname1=[Directory_of_TracesFiles '/' newfolder  '/hel' FileIndexNumber ' tr' num2str(TracesCounter) '-' num2str(Cut_Counter) '.dat'];
            Output_to_be_saved_inaFile=[TimeSeries(x(1):x(2))' Donors(TracesCounter,x(1):x(2))' (Acceptors(TracesCounter,x(1):x(2))-ChannelLeakage*Donors(TracesCounter,x(1):x(2)))'];
            save(fname1,'Output_to_be_saved_inaFile','-ascii') ;
            Done_Cutting_Choice=input('\nPress 0 to continue cutting the same trace \nPress Enter to move to next trace\n');
        end
    end
    
    if choice==4
        TheFileThatWillbeAnalysed_truncated=TheFileThatWillbeAnalysed(1:end-7);
        FilenameForTheImage_Saving=sprintf('%s_Trace_%d.png',TheFileThatWillbeAnalysed_truncated,TracesCounter);
        print(FilenameForTheImage_Saving,'-dpng','-r500');
    end
    
    if choice==5
        break
    end
    
    if choice==6
        while true
            disp('    Click for beginning and end of states.');disp('    Left/middle/right click for different states.');
            [time,y,button]=ginput;
                        
            seld=size(time);
            hold on
            plot(time, y, 'x', 'Color', 'b');
            temp_choice = input(sprintf('You selected %d points\nEnter to accept and go to next image\nPress 9 to re-select\n', seld(1)));
            if temp_choice == 9
                plot(time, y, 'x', 'Color', 'w');
                continue
            end
            
            if mod(seld(1),2)
                disp('You clicked an odd number of times! Please select again');
                plot(time, y, 'x', 'Color', 'w');
                continue;
            end
            
            hold off
            
            % modified to analyze photobleaching rate
            
            %t1 = 10*Timeunit; % 10 red at the beginning
%             t1 = 0;
%             t2 = time(1);
%             DT1(end+1) = abs(t2-t1);
%             fprintf('Dwell time: %f\n', DT1(end));
            
            % original dwell time code
            time1=time(button==1);
            
            for c=1:2:(sum(button==1)-1)
                t1=ceil(time1(c)/Timeunit);
                t2=ceil(time1(c+1)/Timeunit);
                dt1 = abs(time1(c+1) - time1(c));
                fprintf('dt1 = %.2f s\n', dt1);
                DT1(end+1) = dt1;
                DT1a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
                DT1d(end+1)=mean(Donors(TracesCounter,t1:t2));
                DT1f(end+1)=mean(FRET_Time_Series(t1:t2));
            end
            time2=time(button==2);
            for c=1:2:sum(button==2)-1
                t1=ceil(time2(c)/Timeunit);
                t2=ceil(time2(c+1)/Timeunit);
                dt2 = abs(time2(c+1) - time2(c));
                fprintf('dt2 = %.2f s\n', dt2);
                DT2(end+1) = dt2;
                DT2a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
                DT2d(end+1)=mean(Donors(TracesCounter,t1:t2));
                DT2f(end+1)=mean(FRET_Time_Series(t1:t2));
            end
            
            time3=time(button==3);
            for c=1:2:sum(button==3)-1
                t1=ceil(time3(c)/Timeunit);
                t2=ceil(time3(c+1)/Timeunit);
                dt3 = abs(time3(c+1) - time3(c));
                fprintf('dt3 = %.2f s\n', dt3);
                DT3a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
                DT3d(end+1)=mean(Donors(TracesCounter,t1:t2));
                DT3f(end+1)=mean(FRET_Time_Series(t1:t2));
            end
            break
        end
    end

    if choice==7
        choice=input('Which # trace do you want to go to?\n');
        TracesCounter=choice-1;
        continue;
    end
    
%     if choice == 8
%         step_count = input('Number of photobleaching steps in this trace:\n');
%         if ~isempty(step_count)
%             step_count_list(TracesCounter) = step_count;
%         end
%     end

%     if choice == 8
%         [x,y]=ginput(2);
%         st = floor((x(1)-TimeSeries(1))/Timeunit);
%         en = floor((x(2)-TimeSeries(1))/Timeunit);
%         sele_FRET = FRET_Time_Series(st:en);
%         figure;
%         autocorr(sele_FRET);
%         TracesCounter = TracesCounter - 1;
%     end
  
    if choice == 8
        [Xc, ~, ~] = ginput(2);
        start_pt = round(Xc(1)/Timeunit);
        end_pt = round(Xc(2)/Timeunit);
        disp(start_pt);
        disp(end_pt);
        TracesCounter = TracesCounter - 1;
        continue;
    end
    %     if choice == 9
    %         counter = counter + 1;
    %         fprintf('Current counter = %d\n', counter);
    %     end
    
    if choice == 9
        
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
        ylim([-0.2 1.2]);
        
        done_smooth = input('Enter to move on to next trace\n');
      
    end
    
    start_pt = 1;
    end_pt = Length_of_the_TimeTraces;
    
end

%fprintf('counter = %d\n', counter);

% analyze photobleaching step counts
instr = input('Do you want to analyze photobleaching steps? (Y/N, default Y)\n', 's');

if isempty(instr)
    instr = 'Y';
end

if strcmp(instr, 'Y') == 1
    
    figure;
    % create a histogram for photobleaching steps
    count_hist = histogram(step_count_list);
    ch_values = count_hist.Values();
    max_count = size(ch_values, 2);
    xlabel('photobleaching steps');
    ylabel('number of spots');
    
    fname = strcat('hel', FileIndexNumber, '_photobleaching_steps');
    png_fname = strcat(fname,'.png');
    saveas(count_hist, png_fname);
    
    % print the histogram results
    for i_count = 1: max_count
        fprintf('number of %d-step photobleaching: %d\n', ...
            i_count, ch_values(i_count));
    end
    
    % save the counts as a .csv file
    indices = 1:1:num_molecules;
    print_count_list = [indices' step_count_list];
    csv_fname = strcat(fname, '.csv');
    csvwrite(csv_fname, print_count_list);
    
    %     fprintf(fileID,'%10s %10s\n','molecule','step count');
    %     fprintf(fileID, '%10d %10d\n',print_count_list);
    
end

fprintf('Saving dwell time data if there is any...\n');

if ~isempty(DT1)
    DT1=[DT1;DT1a;DT1d;DT1f]';
    fname1=[Directory_of_TracesFiles '/' newfolder  '/dwelltime1.dat'];
    save(fname1,'DT1','-ascii','-append');
end
if ~isempty(DT2)
    DT2=[DT2;DT2a;DT2d;DT2f]';
    fname1=[Directory_of_TracesFiles '/' newfolder  '/dwelltime2.dat'];
    save(fname1,'DT2','-ascii','-append');
end
if ~isempty(DT3)
    DT3=[DT3;DT3a;DT3d;DT3f]';
    fname1=[Directory_of_TracesFiles '/' newfolder  '/dwelltime3.dat'];
    save(fname1,'DT3','-ascii','-append');
end

fprintf('Done.\n');

%close all;
%clear all;

% This will zip all the traces that you wished to convert into pictures.
% Prevents the cluttering of too many files in a folder
% zip('Selected Traces in Images','*.png');
% delete('*.png');


% This will zip all the traces that you selected and saved in the .DAT file format.
% Prevents the cluttering of too many files in a folder
% zip('Selected Traces in DAT files','*.dat');
% delete('*.dat');


