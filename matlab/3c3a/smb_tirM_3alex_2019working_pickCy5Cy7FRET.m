function smb_tirM_3color_2022_pickd2a_XF
% Single Molecule Biophysics Lab. in Seoul National University // MJ 2019 July
% last edited by X. Feng Jan 13, 2022

clear all;
close all;

%% path and filename setting
WorkingDirectory=pwd
filename_head = 'hel1'

%% Correction parameters for FRET%%
dbackground_0=0;
d2background_0=0;
abackground_0=0;

dbackground_1=0;
d2background_1=0;
abackground_1=0;

dbackground_2=0;
d2background_2=0;
abackground_2=0;

leakage01=0.1066;   %0.11
leakage10=0.0;
leakage02=0.0083;   %0.013
leakage12=0.0446;   %0.12

gamma01=0.8730;  %1
gamma12 = 2.62;
gamma02=gamma01*gamma12;

direct = 0.12;

%% Options
LaserOrderChange = 'y'; %Check this part when excitation laser order is matched.
ColorNumber = 3;
DyeType = 'cy235';      %not done: b->g, g->r, r->far red (color change only)
DoseBinningNeed = 'n';  %binning required??
binwidth=5;
DoesFilterNeed = 'n';   %filter required??
DoseMovieNeed = 'n';    % Movie required??
Is_Avg_and_E_save ='n'; % Average and E level save??
Time_unit_ms = 'n';     % Time unit is 'ms'?
FirstNumber = 50;       % histogram options
LastNumber = 10;        % histogram options

%% Trace axis range
BottomLimit_b=-100;
UpperLimit_b=1200;
BottomLimit_g=-100;
UpperLimit_g=1500;
BottomLimit_r=-100;
UpperLimit_r=1500;

%% log file loading (time unit etc. Ha lab ver.)

fileinfo = dir([filename_head '.log']);
if sum(size(fileinfo)) == 1
    disp(['No log file : '  filename_head '.log']);
end
%date = fileinfo(1).date;
fileid_log = fopen([filename_head '.log'],'r');		%% .log file
logarray = textscan(fileid_log, '%s', 'Delimiter','\n');
timeunit = 0.001*str2double(logarray{1,1}{strmatch('Exposure Time [ms]', logarray{1,1})+1})
gain = str2double(logarray{1,1}{strmatch('Gain', logarray{1,1})+1})
scaler = str2double(logarray{1,1}{strmatch('Data Scaler', logarray{1,1})+1})
background_donor = str2double(logarray{1,1}{strmatch('Background', logarray{1,1})+1})
background_acceptor = str2double(logarray{1,1}{strmatch('Background', logarray{1,1})+1})
fclose(fileid_log);

%% path and filename

OriginalDirectory=cd;
cd(WorkingDirectory);
filename_traces = [filename_head '.' num2str(ColorNumber) 'color_3alex_traces'];
filename_movie = [filename_head '.' num2str(ColorNumber) 'color_3alex_movies'];
filename_time = [filename_head '_time.dat'];
filename_add = [ filename_head '_select.dat'];
filename_add_region = [ filename_head '_region.dat' ];
filename_add_region_data=[ filename_head '_region_data.dat' ];

%% Reading data %%
fileid = fopen(filename_traces,'r');

if fileid == -1
    disp(['No data file  '  filename_head]);
end

time_length = fread(fileid, 1, 'int32');
disp(sprintf('The length of the time traces is: %d', time_length));

time = (0:(time_length-1))*timeunit;

NumberofTraces = fread(fileid, 1, 'int16');
NumberofPeaks = NumberofTraces/ColorNumber;
disp('The number of traces and peaks are:')
disp(NumberofTraces);
disp(NumberofPeaks);
%Data = zeros(NumberofTraces, time_length, 'int16');
Data = fread(fileid, [NumberofTraces  time_length],'int16');
SpotDiameter = fread(fileid, 1, 'int16');
disp('Done reading trace data.');
fclose(fileid);

% Temp_i = [];
% Tempfirstpoint = [];
% Templastpoint = [];
% Templength = [];
% Temptime_region = [];
% TempFret12_region = [];
% TempFret13_region = [];
% TempFret23_region = [];
% TempDonor_b_region = [];
% TempDonor2_b_region = [];
% TempAcceptor_b_region = [];
% TempDonor_g_region = [];
% TempDonor2_g_region = [];
% TempAcceptor_g_region = [];
% TempDonor_r_region = [];
% TempDonor2_r_region = [];
% TempAcceptor_r_region = [];
dbackground_0_temp = dbackground_0;
d2background_0_temp = d2background_0;
abackground_0_temp = abackground_0;

fileid_add_region = fopen(filename_add_region, 'r');
%fileid_add_region_data = fopen(filename_add_region_data, 'r');
fileid_add_region_data = -1;
if fileid_add_region~=-1 && fileid_add_region_data~=-1;
    fgetl(fileid_add_region);
    temp = fscanf(fileid_add_region, '%f %f %f  %f %f %f  %f %f %f  %f %f %f %f  %f %f  %f\n', [13 1]);
    
    dbackground_0=double(temp(1,1));
    d2background_0=double(temp(2,1));
    abackground_0=double(temp(3,1));
   
    
    dbackground_1=double(temp(4,1));
    d2background_1=double(temp(5,1));
    abackground_1=double(temp(6,1));
    
    dbackground_2=double(temp(7,1));
    d2background_2=double(temp(8,1));
    abackground_2=double(temp(9,1));
    
    leakage01=double(temp(10,1));   %0.11
    leakage10=double(temp(11,1));
    leakage02=double(temp(12,1));   %0.013
    leakage12=double(temp(13,1));   %0.12
    
    gamma01=double(temp(14,1));  %1
    gamma02=double(temp(15,1));  %2.36
    
    direct=double(temp(16,1));   %0.19
    
    fgetl(fileid_add_region);
    temp = fscanf(fileid_add_region, '%f %f %f %f\n', [4 inf]);
    Temp_i = int32(temp(1,:));
    Tempfirstpoint = int32(temp(2,:));
    Templastpoint = int32(temp(3,:));
    Templength = int32(temp(4,:));
    
    temp = fscanf(fileid_add_region_data, '%f %f %f %f  %f %f %f  %f %f %f  %f %f %f\n', [10 inf]);
    Temptime_region = temp(1,:);
    TempFret12_region = temp(2,:);
    TempFret13_region = temp(3,:);
    TempFret23_region = temp(4,:);
    TempDonor_b_region = temp(5,:);
    TempDonor2_b_region = temp(6,:);
    TempAcceptor_b_region = temp(7,:);
    TempDonor_g_region = temp(8,:);
    TempDonor2_g_region = temp(9,:);
    TempAcceptor_g_region = temp(10,:);
    TempDonor_r_region = temp(11,:);
    TempDonor2_r_region = temp(12,:);
    TempAcceptor_r_region = temp(13,:);
    fclose(fileid_add_region);
    fclose(fileid_add_region_data);
end


if DoseMovieNeed == 'y'
    cd(WorkingDirectory);
    disp(filename_movie);
    fileid_movie = fopen(filename_movie, 'r');
    if fileid_movie ~= -1
        peaks_total_width = fread(fileid_movie, 1, 'int16');
        peak_height = fread(fileid_movie, 1, 'int16');
        peaks_number = peaks_total_width/peak_height;
        file_information = dir(filename_movie);
        film_time_length = (file_information.bytes-4)/(peak_height*peaks_total_width);
        fclose(fileid_movie);
        
        disp('peaks_total_width, height, number, film_time_length: ');
        disp(peaks_total_width);
        disp(peak_height);
        disp(peaks_number);
        disp(film_time_length);
        
        peak_line=zeros(ColorNumber*peak_height, 1, 'uint8');
        %%peaks=fread(fileid_movie, peaks_total_width * peak_height * film_time_length, 'uint8');
        peak=zeros(ColorNumber*peak_height, peak_height, 'uint8');
        
        if peaks_number ~= NumberofTraces
            disp('error: Different trace numbers between .trace and .movies');
            return;
        end
        
        if film_time_length ~= time_length
            disp('error: Different time length between .trace and .movies');
            return;
        end
    else
        DoseMovieNeed = 'n';
        disp('No movie file');
    end
end



%% Convert raw data into donor and acceptor traces %%
time_length_each = time_length;

DonorRawData = zeros(NumberofPeaks, time_length_each, 'double');
Donor2RawData = zeros(NumberofPeaks, time_length_each, 'double');
AcceptorRawData = zeros(NumberofPeaks, time_length_each, 'double');

binlength = int32(time_length/binwidth-1);
bintime = zeros(binlength, 1, 'double');
binEraw = zeros(binlength, 1, 'double');
binEcorrect = zeros(binlength, 1, 'double');


for m=1:binlength
    bintime(m) = double(m-1)*(binwidth*timeunit);
end



for i = 1:NumberofPeaks
    
    DonorRawData(i, :) = Data(i*3-2, :);
    Donor2RawData(i, :) = Data(i*3-1, :);
    AcceptorRawData(i, :) = Data(i*3, :);
    
end


clear Data;



%% calculate, plot and save average traces %%
%MJ edited

DonorRawData_0 = DonorRawData;
Donor2RawData_0 = Donor2RawData;
AcceptorRawData_0 = AcceptorRawData;

DonorRawData_1 = DonorRawData(:, 13:18);
Donor2RawData_1 = Donor2RawData(:, 13:18);
AcceptorRawData_1 = AcceptorRawData(:, 13:18);

%% Start to servey
Temptime = [];
%TempFret = [];
%TempDonor = [];
%TempAcceptor = [];
TempFret12 = [];
TempFret13 = [];
TempFret23 = [];
TempDonor_b = [];
TempDonor2_b = [];
TempAcceptor_b = [];
TempDonor_g = [];
TempDonor2_g = [];
TempAcceptor_g = [];
TempDonor_r = [];
TempDonor2_r= [];
TempAcceptor_r = [];
r12_list = [];
r13_list = [];
r23_list = [];
l12_list = [];
l13_list = [];
l23_list = [];
dd2ac_list = [];
t_list = -1 * ones(NumberofPeaks, 1);
DT1=[];DT2=[];DT3=[];
DT1start = [];
DT1end = [];
Cy57_avg_FRET_list = zeros(NumberofPeaks, 1);
% DT1a=[];DT2a=[];DT3a=[];
% DT1d=[];DT2d=[];DT3d=[];
% DT1f=[];DT2f=[];DT3f=[];
DT1id = [];
dtID = 0;
dtCell = cell(0, 6);

scrsz = get(0,'ScreenSize');
figure('Name','Trace analysis','OuterPosition',[scrsz(4)*0.1 0.1*scrsz(4) scrsz(4)*1.5 0.8*scrsz(4)]);

i=0;
history_n = 0;
history = zeros(1000, 1, 'int16');

newfolder = [filename_head '_analysis'];
mkdir(newfolder);

% Display traces

while i < NumberofPeaks
    disp(i);
    i = int32(int32(i) + 1);
    
    again=1;
    while ((again==1) && (i <= NumberofPeaks))
        if ColorNumber == 3
            if DyeType == 'cy235'
                DonorCorrect_0 = ((1 + leakage01 + leakage02) / (1 - leakage01 * leakage10)) * ((DonorRawData_0(i,:) - dbackground_0_temp) - leakage10 * (DonorRawData_0(i,:) - d2background_0_temp));
                Donor2Correct_0 = gamma01 * ((1 + leakage10 + leakage12) / (1 - leakage01 * leakage10)) * ((Donor2RawData_0(i,:) - d2background_0_temp) - leakage01 * (DonorRawData_0(i,:) - dbackground_0_temp));
                AcceptorCorrect_0 = gamma02 * ((AcceptorRawData_0(i,:) - abackground_0_temp) - ((leakage12 * leakage01 - leakage02)/(1 - leakage01 * leakage10)) * (DonorRawData_0(i,:) - dbackground_0_temp) + ((leakage02 * leakage10 - leakage12 ) / (1 - leakage01 * leakage10)) * (Donor2RawData_0(i,:) -d2background_0_temp));
                EachTotalCorrect_0 = DonorCorrect_0 + Donor2Correct_0 + AcceptorCorrect_0;
                EachTotalCorrect_0 = (EachTotalCorrect_0~=0).*EachTotalCorrect_0 + (EachTotalCorrect_0==0)*1;	% remove zeros
                
                %DonorCorrect_g = ((1 + leakage01 + leakage02) / (1 - leakage01 * leakage10)) * ((DonorRawData_1(i,:) - dbackground_1) - leakage10 * (Donor2RawData_1(i,:) - d2background_1));
                Donor2Correct_1 = gamma01 * ((1 + leakage10 + leakage12) / (1 - leakage01 * leakage10)) * ((Donor2RawData_1(i,:) - d2background_1) - leakage01 * (DonorRawData_1(i,:) - dbackground_1));
                AcceptorCorrect_1 = gamma02 * ((AcceptorRawData_1(i,:) - abackground_1) - ((leakage12 * leakage01 - leakage02)/(1 - leakage01 * leakage10)) * (DonorRawData_1(i,:) - dbackground_1) + ((leakage02 * leakage10 - leakage12 ) / (1 - leakage01 * leakage10)) * (Donor2RawData_1(i,:) -d2background_1));
                AcceptorCorrect_1 = AcceptorCorrect_1 - direct * (Donor2Correct_1 + AcceptorCorrect_1);
                EachTotalCorrect_1 = Donor2Correct_1 + AcceptorCorrect_1;
                EachTotalCorrect_1 = (EachTotalCorrect_1~=0).*EachTotalCorrect_1 + (EachTotalCorrect_1==0)*1;	% remove zeros
                
                Fret12 = mean(AcceptorCorrect_1./EachTotalCorrect_1);
                %disp(AcceptorCorrect_1);
                %disp(EachTotalCorrect_1);
                Fret01 = Donor2Correct_0./((1-Fret12)*DonorCorrect_0 + Donor2Correct_0);
                Fret02 = (AcceptorCorrect_0 - Fret12*(Donor2Correct_0 + AcceptorCorrect_0))./(DonorCorrect_0 + AcceptorCorrect_0 - Fret12 * (EachTotalCorrect_0));
            
                Cy57_avg_FRET_list(i) = Fret12;
                %disp(Cy57_avg_FRET_list(i));
            end
        end
        
        
        
        for j=2:time_length_each
            if (Fret02(j) < -0.3)
                Fret02(j) = nan;
            elseif (Fret02(j) > 1.3)
                Fret02(j) = nan;
            end
            
            if(Fret01(j) < -0.3)
                Fret01(j) = nan;
            elseif(Fret01(j) > 1.3)
                Fret01(j) = nan;
            end
%             if(Fret12(j) < -0.3)
%                 Fret12(j) = 0;
%             elseif(Fret12(j) > 1.1)
%                 Fret12(j) = 1;
%             end

            if EachTotalCorrect_0(j) < 100
                Fret02(j) = nan;
                Fret01(j) = nan;
            end
        end
        
        
        
        again=0;
        %end
        %i = i + 1 ;
    end
    
    subplot(3, 1, 1);
    plot(time, DonorCorrect_0, 'Color', [0.13, 0.53, 0.13]);
    hold on
    plot(time, Donor2Correct_0, 'r', time, AcceptorCorrect_0, 'b');
    hold off
%     hold on
    %plot(time, DonorCorrect_0 + Donor2Correct_0 + AcceptorCorrect_0 + 300, 'k');
%     hold off
    temp=axis;
    temp(3)=BottomLimit_b;
    temp(4)=UpperLimit_b;
    grid on;
    axis(temp);
    %title(['                    Molecule ' num2str(i) '  / ' num2str(NumberofPeaks) '              file: ' filename_head ]); %'   cor.: ' num2str(dbackground) '  ' num2str(abackground) '  ' num2str(d2background) ' ' num2str(leakage01) '  ' num2str(leakage02) ' ' num2str(leakage12) ' ' num2str(gamma01) ' ' num2str(gamma02) '   filter: ' DoesFilterNeed ' Fret:' num2str(FretEc_filter) ' Mean:' num2str(MeanRange) ' Acceptor:' num2str(Acceptor_filter) ' Total:' num2str(Total_filter) ]);
    title(['ALEX Green Laser Molecule' num2str(i) '  / ' num2str(NumberofPeaks) '       file:' filename_head  ' cor.: ' num2str(dbackground_0) '  ' num2str(d2background_0) '  ' num2str(abackground_0) '  ' num2str(dbackground_1) '  ' num2str(d2background_1) '  ' num2str(abackground_1) '  ' num2str(dbackground_2) '  ' num2str(d2background_2) '  ' num2str(abackground_2) '  ' num2str(leakage01) '  ' num2str(leakage10) '  ' num2str(leakage02) '  ' num2str(leakage12) '  ' num2str(gamma01) '  ' num2str(gamma02) ]);
%     zoom on;
    
    subplot(3, 1, 2);
    plot(time, Fret01, 'r');
    ylim([-0.2 1.2]);
    title('Cy3 Cy5 FRET');
    
    subplot(3, 1, 3);
    plot(time, Fret02, 'b');
    ylim([-0.2 1.2]);
    title('Cy3 Cy7 FRET');
    
%     disp([num2str(i) ' (d=collect dwell times)']);
    disp('s = save subtrace, d = collect dwelltime');
    keyanswer =input('(t=terminate program, b=back, g=go, o=subtract background, q=collect photobleaching time) : ','s');
    answer = sscanf(keyanswer, '%s %*s');
    numberofanswer = sscanf(keyanswer, '%*s %f');
%     answer = 0;

    % photobleaching analysis
    if answer == 'q'
        [x, y] = ginput(1);
        disp(x);
        t_list(i) = x;
    end
    
    % dwell time analysis
    if answer == 'd'
        disp('Click for beginning and end of states.');disp('Left/middle/right click for different states.');
        [times,y,button]=ginput;
        
        time1=times(button==1);
        for c=1:2:(sum(button==1)-1)
            
            t1=ceil(time1(c)/timeunit);
            t2=ceil(time1(c+1)/timeunit);
            DT1(end+1)=abs(time1(c+1)-time1(c));
            DT1id(end+1) = i;
            DT1start(end+1) = time1(c);
            DT1end(end+1) = time1(c+1);
            
            dtID = dtID + 1;
            dtCell{dtID, 1} = i;
            dtCell{dtID, 2} = time1(c);
            dtCell{dtID, 3} = time1(c+1);
            dtCell{dtID, 4} = Fret12;
            dtCell{dtID, 5} = Fret01(t1:t2);
            dtCell{dtID, 6} = Fret02(t1:t2);
            %             DT1a(end+1)=mean(Acceptors(i,t1:t2));
            %             DT1d(end+1)=mean(Donors(i,t1:t2));
            %             DT1f(end+1)=mean(FRET_Time_Series(t1:t2));
        end
        time2=times(button==2);
        for c=1:2:sum(button==2)-1
            t1=ceil(time2(c)/timeunit);t2=ceil(time2(c+1)/timeunit);
            DT2(end+1)=abs(time2(c+1)-time2(c));
            %             DT2a(end+1)=mean(Acceptors(i,t1:t2));
            %             DT2d(end+1)=mean(Donors(i,t1:t2));
            %             DT2f(end+1)=mean(FRET_Time_Series(t1:t2));
        end
        time3=times(button==3);
        for c=1:2:sum(button==3)-1
            t1=ceil(time3(c)/timeunit);t2=ceil(time3(c+1)/timeunit);
            DT3(end+1)=abs(time3(c+1)-time3(c));
            %             DT3a(end+1)=mean(Acceptors(i,t1:t2));
            %             DT3d(end+1)=mean(Donors(i,t1:t2));
            %             DT3f(end+1)=mean(FRET_Time_Series(t1:t2));
        end
        
    end
    
    if answer == 'o'
            
            [raw_x, ~] = ginput(2);
            x = round(raw_x / timeunit);
            dbackground_0_temp = mean(DonorRawData_0(i,x(1):x(2)));
            d2background_0_temp = mean(Donor2RawData_0(i, x(1):x(2)));
            abackground_0_temp = mean(AcceptorRawData_0(i, x(1):x(2)));
            
            % std
            dbackground_0_std = std(DonorRawData_0(i,x(1):x(2)));
            d2background_0_std = std(Donor2RawData_0(i, x(1):x(2)));
            abackground_0_std = std(AcceptorRawData_0(i, x(1):x(2)));
            
            output = [dbackground_0_temp; d2background_0_temp; abackground_0_temp; dbackground_0_std; d2background_0_std; abackground_0_std];
            file_outname = sprintf('%s/%s_trace_%i_background.dat', newfolder, filename_head, i);
            save(file_outname,'output','-ascii');
            
            i = i - 1;
            
            continue;
    end
    
    if answer=='b'
        if i>1 && history_n > 0
            while history(history_n)==i
                history_n = history_n - 1;
            end
            i=history(history_n)-1;
            history_n = history_n - 1;
        else
            i=i-1;
        end
    else
        history_n = history_n + 1;
        history(history_n)=i;
    end
    
    
    if answer=='g'
        answer=input('number to go : ','s');
        gonumber = str2num(answer);
        if gonumber > 0 && gonumber <= NumberofPeaks
            i = gonumber - 1;
        end
    end
    
    if answer=='t'
        close all;
        cd(OriginalDirectory);
        break;
    end
    
    if answer=='s'
        [Xc,~,buttonc] = ginput(2); 
        if (Xc(1)>Xc(2))
            temp = Xc(1);
            Xc(1) = Xc(2);
            Xc(2) = temp;
        end
        firstpoint = round(Xc(1)/timeunit); 
        lastpoint = round(Xc(2)/timeunit); 
        disp(firstpoint);
        disp(lastpoint);
        Temptime = time(firstpoint:lastpoint);
        TempFret01 = Fret01(firstpoint:lastpoint);
        TempFret02 =  Fret02(firstpoint:lastpoint);
        Temp0 = DonorCorrect_0(firstpoint:lastpoint);
        Temp1 = Donor2Correct_0(firstpoint:lastpoint);
        Temp2 = AcceptorCorrect_0(firstpoint:lastpoint);
       
        output = [Temptime' Temp0' Temp1' Temp2' TempFret01' TempFret02'];
        file_outname = sprintf('%s/%s_trace_%i.dat', newfolder, filename_head, i);
        save(file_outname,'output','-ascii');
    end
    
    dbackground_0_temp = dbackground_0;
    d2background_0_temp = d2background_0;
    abackground_0_temp = abackground_0;
    
end

disp('Program end.')

cd(OriginalDirectory);

% save dwell time data

fprintf('Saving dwell time data if there is any...\n');


% if ~isempty(DT1)
%     %DT1=[DT1;DT1a;DT1d;DT1f]';
%     DT1table = [DT1id; DT1start; DT1end; DT1]';
%     fname1=[WorkingDirectory '/' newfolder  '/dwelltime1.dat'];
%     save(fname1,'DT1table','-ascii','-append');
% end
% if ~isempty(DT2)
%     DT2=[DT2;DT2a;DT2d;DT2f]';
%     fname1=[Directory_of_TracesFiles '/' newfolder  '/dwelltime2.dat'];
%     save(fname1,'DT2','-ascii','-append');
% end
% if ~isempty(DT3)
%     DT3=[DT3;DT3a;DT3d;DT3f]';
%     fname1=[Directory_of_TracesFiles '/' newfolder  '/dwelltime3.dat'];
%     save(fname1,'DT3','-ascii','-append');
% end

% TraceIDs = 1:NumberofPeaks;
% Cy57_fret_table = [TraceIDs' Cy57_avg_FRET_list];
% fname_57fret = [WorkingDirectory '/' newfolder  '/cy57_fret.dat'];
% save(fname_57fret, 'Cy57_fret_table', '-ascii');

if size(dtCell, 1) > 0
    fname = [newfolder '/selected_raw_FRET.mat'];
    col_names = ["TraceID" "Start" "End" "Cy3Em" "Cy5Em" "Cy7Em" "Fret12" "Fret01" "Fret02"];
    save(fname, 'dtCell', 'col_names');
end

fprintf('Done.\n');

clear all;
close all;
end
