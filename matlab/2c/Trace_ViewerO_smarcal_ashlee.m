% Written by Digvijay Singh  ( dgvjay@illinois.edu)
% This script will let you the single molecule traces from the .traces file
% After entering the required information initially.
% It will sift through each trace present in the analysed movie file ( .traces file)
% and you can decide to SAVE THEM in different ways...read on to
% understand or simply RUN this script and you will know everything
% Edited by Olivia Yang (oyang1@jhu.edu) and Ashlee Feng (xfeng17@jhu.edu)

%close all;
clear
fclose('all');
minInt=100;
ymax=800;
Directory_of_TracesFiles=input('Directory: ','s');

if isempty(Directory_of_TracesFiles)
    Directory_of_TracesFiles=pwd;
end

cd(Directory_of_TracesFiles);

% This script ases that the .traces file in the folder(input above) is
% in the following order in your directory
% like hel1.traces...hel2.traces...hel3.traces

GenericFileType='.traces';   % The single molecule information post IDL analysis is stored in the .traces file
% for each movie analysed.

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

% We typically employ a GammaFactor of 1.0 but if you want to change
% it...you can change it with the routine below
GammaFactor=input('Enter the value of the Gamma Factor [Default is 1.0] :  ');

if isempty(GammaFactor)
    GammaFactor=1.0;
end

% Part of the FRET Calculation
% Channel Leakage Value
ChannelLeakage=input('Enter the value of the Channel Leakage [Default is 0.12] ');
if isempty(ChannelLeakage)
    ChannelLeakage=0.12;
    % T70S: 0.175
end


% Extracting important information from .traces binary file
Length_of_the_TimeTraces=fread(File_id,1,'int32');% This fetches the total duration of the imaging that was carried out for the

% concerned .traces file...please note that .traces is a binary file and
% this is the way it was binarized and we are just extracting the
% information from the binary file
disp('The length of the time traces is: ')
disp(Length_of_the_TimeTraces);% This displays the total duration of the imaging that was carried out for the
% concerned .traces file.

num_traces=fread(File_id,1,'int16');  % This fetches the total number of single molecule spots in the
% concerned .traces file...since each spot has a pair i.e. Green channel
% and a red channel therefore we would need to divide this number by 2 to
% get the actual number of spots.

disp('The number of traces in this file is:')
num_molecules = num_traces / 2;
disp(num_molecules);% This displays the total number of single molecule spots in the
% concerned .traces file.

%Reading in the entire raw data from the .traces binary file encoded in
%int16 format as noted above that .traces is a binary file and
% this is the way it was binarized and we are just extracting the
% information from the binary file.
Raw_Data=fread(File_id,num_traces*Length_of_the_TimeTraces,'int16');
disp('Done reading data');
fclose(File_id);  % We close the file pointer here as all the information that we needed from it
% has been extracted succesffuly and stored into local variables like Raw_Data, num_traces etc etc.


% Converting into Donor and Acceptor traces of several selected spots in
% the movie.
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
%TimeSeries is nothing but the index of each frame event in the series.
% Suppose you made a movie which is 10 seconds long, if the frame rate was 100ms. It means you have collected
%     10/0.1 ==100 frames. TimeSeries is simply a series from 1 to 100 where index 1 would correspond to first frame
%     index 2 would correspond to frame2 and index 3 would correspond to frame numbrr 3 and so on and each frame number
%     would have a particular intensity for the spot in the Green channel and the red channel.

smoothed_fret_x = zeros(floor(Length_of_the_TimeTraces/3));
for i = 1: floor(Length_of_the_TimeTraces/3)
    smoothed_fret_x(i) = TimeSeries(i*3);
end

% open and read the .pks for coordinates
% ['hel' FileIndexNumber '.pks'] - unnecessary?
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


%select the folder to which the files need to be saved
newfolder = [FileIndexNumber ' selected traces'];
mkdir(newfolder);

% check for constructs with both dyes
%HasBoth = sum(mean(Acceptors(:,end-10:end-1)')>=400)

% ==========================================================================

TracesCounter=0;
%Now we will be going over the trace one by one to select them or whatever
figure;
hd13=gcf;

%flowtime=100*Timeunit;

DT1=[];DT2=[];DT3=[];
DT1a=[];DT2a=[];DT3a=[];
DT1d=[];DT2d=[];DT3d=[];
DT1f=[];DT2f=[];DT3f=[];

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
    ax1=subplot(2,4,[1 2 3]);
    %     ax1 = subplot(2, 1, 1);
    %Simply plotting the Donor and the Leakage Corrected Acceptor
    %    plot(TimeSeries,Donors(TracesCounter,:),'g',...
    %        TimeSeries,(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:)),'r',...
    %        TimeSeries,300+Donors(TracesCounter,:)+(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:)),'k');
    
    plot(TimeSeries,Donors(TracesCounter,:),'g',...
        TimeSeries,(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:)),'r');
    
    %FLOW
    %line([flowtime,flowtime],[-5000,5000])
    
    %HIGH FRET  only
    %if mean(Acceptors(TracesCounter,1:20))<400
    %    continue;
    %end
    
    %HIGH Cy5 end frames  only
    %if mean(Acceptors(TracesCounter,end-10:end-1))<400
    %    continue;
    %end
    
    %HIGH Cy5 beginning frames  only
    %if mean(Acceptors(TracesCounter,1:10))<400
    %    continue;
    %end
    
    %hold
    %plot(TimeSeries,(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:)),'LineWidth',1.2,'Color','r');
    %Turn the below on if you also want to see the total intensity time series( i.e.
    %Acceptor Plus Donor)
    %legend({'Donor Intensity ','Acceptor Intensity'},'FontSize',8,'FontWeight','bold')
    temp=axis;
    temp(2)=TimeSeries(end);
    temp(3)=-100;
    %    temp(4)=ymax; %adjust max y-axis
    axis(temp);
    
    
    xlabel('Time(s)');
    ylabel('Intensity (a.u.)');
    xlim([0 Length_of_the_TimeTraces * Timeunit]);
    %plot(TimeSeries,Acceptors(TracesCounter,:)+Donors(TracesCounter,:),'k');
    TitleNameForThePlot=sprintf('Molecule %d / %d of %s',TracesCounter,num_traces/2,TheFileThatWillbeAnalysed);
    title(TitleNameForThePlot); % Giving the title name to the plot as described above
    
    ax2=subplot(2, 4,[5 6 7]);
    %     ax2 = subplot(2, 1, 2);
    %Simply plotting the  FRET traces now in a subplot below the above plot.
    FRET_Time_Series=(Acceptors(TracesCounter,:)...
        -ChannelLeakage*Donors(TracesCounter,:))...
        ./(Acceptors(TracesCounter,:)...
        -ChannelLeakage*Donors(TracesCounter,:)...
        +(Donors(TracesCounter,:)));
    %TitleNameForThePlot=sprintf('Molecule Number %d of the file %s',i,TheFileThatWillbeAnalysed);
    
    %    for m=1:Length_of_the_TimeTraces
    %        if FRET_Time_Series(m) < - 0.2 || FRET_Time_Series(m) > 1.2
    %            FRET_Time_Series(m)=NaN;
    %        end
    %    end
    
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
    
    %     subplot(3,3,7);
    %     hist(FRET_Time_Series,-.1:.025:1.1);
    %     xlabel('E_F_R_E_T'); ylabel('Count'); title('Trace FRET histogram');
    %     temp=axis;   temp(1)=-0.1;   temp(2)=1.1;   axis(temp);
    %
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
        subplot(2, 4, 8);
        imshow(imresize(acceptorInserter(int16(max(1,aspoty-20)):int16(min(sz,aspoty+20)),int16(max(1,aspotx-20)):int16(min(512,aspotx+20)),:),4,'nearest'));
        title('Acceptor');
        subplot(2, 4, 4);
        imshow(imresize(donorInserter(int16(max(1,dspoty-20)):int16(min(sz,dspoty+20)),int16(max(1,dspotx-20)):int16(min(512,dspotx+20)),:),4,'nearest'));
        title('Donor');
        zoom on;
        
    end
    
    % Now the particular trace has been plotted completely, we now can
    % decide what we want to do with the trace.
    
    % You can choose to
    
    %  PRESS 1.
    %  Save the trace as it is....( At this point it will output a .dat file in the same folder
    %  % where the first column will be Time
    %  the second column will be corrected donor intensitty
    %  the third column will be corrected acceptor intensity
    % (PLEASE NOTE THAT THIS .DAT file's donor and acceptor values will  be corrected for leakage or additional
    %  background.).
    % YOU MAY CHOOSE TO CHANGE THAT BY CHANGING THE CODING HERE.
    
    % PRESS 3.
    % Sometimes you may want to use only a portion of the entire single
    % molecule trace and discard the rest ( for e.g. removing the single
    % step photobleaching event or remove the weird intensity fluctuation
    % regions)
    % In that case you can manually cut out the portion of the trace which is good and relevant and can
    % be used for your studies.
    
    % For cutting out the trace, all you have to do is make two clicks.
    % The first click will be the starting point of where you want to induce the cut and the second click will be the
    % end point
    
    % You may want to make more cuts than just 1 and the script allows you
    % do to that. If you keep pressing 0 after pressing 3...then it will
    % keep letting you make more cuts and all of those cuts will be stored
    % separately in the same format.
    
    % If you want to stop making the cuts, then simply press 1 after making
    % the required number of cuts.
    
    
    % PRESS 0 to go back and look at the previous traces that have gone by
    % and you can choose to press 1 or 3 over them again to select them
    
    
    % PRESS 4 to save the trace in a picture file. You can use this for
    % meeting or your personal records.
    % if you want to save the trace in both .DAT format as well as save
    % image of it...just press 1 or 3 once to save the .DAT file of the
    % trace and then go back to the same trace again by pressing 0 and now
    % press 4 to save it's image....simple !
    
    choice = input('Press 0 to go back one trace\nPress 1 to save current trace\nPress 2 to subtract background (2 clicks)\nPress 3 to cut out\nPress 4 to save the trace as a picture\nPress 5 to skip to data analysis and exit\nPress 6 to collect dwell times\nPress 7 to choose a trace\nPress 8 to count photobleaching steps\nPress 9 to increment counter\nPress Enter to move on\n');
    
    if choice==0
        TracesCounter=TracesCounter - 2;
        continue;
    end
    
    % Making this choice will simply save the entire trace in its format
    % where the first column will be Time
    % the second column will be corrected donor intensitty
    % the third column will be corrected acceptor intensity
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
    
    % Making this choice will simply cut out the trace...when you specify the region
    % The format will be like:
    % where the first column will be Time
    % the second column will be donor intensitty
    % the third column will be raw(non leakage corrected) acceptor intensity
    if choice == 3
        TheFileThatWillbeAnalysed_truncated=TheFileThatWillbeAnalysed(1:end-7);
        fprintf('Click twice to select the range to extract\n');
        Done_Cutting_Choice = 0;
        Cut_Counter = 0;
        while Done_Cutting_Choice == 0
            Cut_Counter = Cut_Counter+1;
            [x, y, button] = ginput(2);  % Make Two clicks to specify a region which you want to cut out.
            x(1) = round(x(1)/Timeunit); % Starting point of the cut region
            x(2) = round(x(2)/Timeunit); % End point of the cut region
            %Only the timeseries and the corresponding Donor and Acceptor intensity falliny betwee the above
            %two points will be stored now, the rest will not be stored.
            %Filename_for_This_Trace=sprintf('%s_Trace_%d_CutCount_%d.dat',TheFileThatWillbeAnalysed_truncated,TracesCounter, Cut_Counter);
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
            % Asking you whether you want to keep cutting the trace or want
            % to move on to the NEXT set of traces.
            Done_Cutting_Choice=input('\nPress 0 to continue cutting the same trace \nPress Enter to move to next trace\n');
        end
    end
    
    % If you press 4....it will save the trace in a picture file. You can use this for
    % meeting or your personal records...
    % if you want to save the trace in both .DAT format as well as save
    % image of it...just press 1 or 3 once to save the .DAT file of the
    % trace and then go back to the same trace again by pressing 0 and now
    % press 4 to save it's image....simple !
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
            
            % for testing
            
            seld=size(time);
            hold on
            plot(time, y, 'x', 'Color', 'b');
            temp_choice = input(sprintf('You selected %d points\nEnter to accept and go to next image\nPress 9 to re-select\n', seld(1)));
            if temp_choice == 9
                plot(time, y, 'x', 'Color', 'w');
                continue
            end
            hold off
            %         if mod(seld(1),2)
            %             disp('Missing one end! Please go back one molecule to select again');
            %             continue;
            %         else
            
            % modified to analyze photobleaching rate
            
            %t1 = 10*Timeunit; % 10 red at the beginning
            t1 = 0;
            t2 = time(1);
            DT1(end+1) = abs(t2-t1);
            fprintf('Dwell time: %f\n', DT1(end));
            
            % original dwell time code
            %                 time1=time(button==1);
            %
            %                 for c=1:2:(sum(button==1)-1)
            %                     t1=ceil(time1(c)/Timeunit);
            %                     t2=ceil(time1(c+1)/Timeunit);
            %                     DT1(end+1)=abs(time1(c+1)-time1(c));
            %                     DT1a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
            %                     DT1d(end+1)=mean(Donors(TracesCounter,t1:t2));
            %                     DT1f(end+1)=mean(FRET_Time_Series(t1:t2));
            %                 end
            %                 time2=time(button==2);
            %                 for c=1:2:sum(button==2)-1
            %                     t1=ceil(time2(c)/Timeunit);t2=ceil(time2(c+1)/Timeunit);
            %                     DT2(end+1)=abs(time2(c+1)-time2(c));
            %                     DT2a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
            %                     DT2d(end+1)=mean(Donors(TracesCounter,t1:t2));
            %                     DT2f(end+1)=mean(FRET_Time_Series(t1:t2));
            %                 end
            %                 time3=time(button==3);
            %                 for c=1:2:sum(button==3)-1
            %                     t1=ceil(time3(c)/Timeunit);t2=ceil(time3(c+1)/Timeunit);
            %                     DT3(end+1)=abs(time3(c+1)-time3(c));
            %                     DT3a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
            %                     DT3d(end+1)=mean(Donors(TracesCounter,t1:t2));
            %                     DT3f(end+1)=mean(FRET_Time_Series(t1:t2));
            %                 end
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

    if choice == 8
        [x,y]=ginput(2);
        st = floor((x(1)-TimeSeries(1))/Timeunit);
        en = floor((x(2)-TimeSeries(1))/Timeunit);
        sele_FRET = FRET_Time_Series(st:en);
        figure;
        autocorr(sele_FRET);
        TracesCounter = TracesCounter - 1;
    end
%     
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
end

fprintf('counter = %d\n', counter);

% analyze photobleaching step counts
instr = input('Do you want to analyze photobleaching steps? (Y/N)\n', 's');

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

close all;
%clear all;

% This will zip all the traces that you wished to convert into pictures.
% Prevents the cluttering of too many files in a folder
% zip('Selected Traces in Images','*.png');
% delete('*.png');


% This will zip all the traces that you selected and saved in the .DAT file format.
% Prevents the cluttering of too many files in a folder
% zip('Selected Traces in DAT files','*.dat');
% delete('*.dat');


