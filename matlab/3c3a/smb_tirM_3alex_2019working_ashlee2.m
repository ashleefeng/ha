function smb_tirM_3alex_ashlee
% Single Molecule Biophysics Lab. in Seoul National University // MJ 2019 July
% edited by X. Feng Sep 4, 2019

clear all;
close all;

%% path and filename setting
WorkingDirectory=pwd
filename_head = 'hel5'

%% Correction parameters for FRET%%
dbackground_b=-50;
d2background_b=0;
abackground_b=0;

dbackground_g=0;
d2background_g=0;
abackground_g=0;

dbackground_r=0;
d2background_r=0;
abackground_r=0;

leakage12=0.1066;   %0.11
leakage21=0.0;
leakage13=0.0083;   %0.013
leakage23=0.0446;   %0.12

gamma12=0.8730;  %1
gamma23 = 2.62;
gamma13=gamma12*gamma23;

direct = 0; %ashlee: 0.1578;

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
UpperLimit_b=1500;
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

if mod(time_length,3)==0
    time_b = (0:+3:(time_length-3))*timeunit;
    time_g = (1:+3:(time_length-2))*timeunit;
    time_r = (2:+3:(time_length-1))*timeunit;
end

if mod(time_length,3)==1
    time_b = (0:+3:(time_length-3))*timeunit;
    time_g = (1:+3:(time_length-5))*timeunit;
    time_r = (2:+3:(time_length-4))*timeunit;
end

if mod(time_length,3)==2
    time_b = (0:+3:(time_length-3))*timeunit;
    time_g = (1:+3:(time_length-2))*timeunit;
    time_r = (2:+3:(time_length-4))*timeunit;
end

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

Temp_i = [];
Tempfirstpoint = [];
Templastpoint = [];
Templength = [];
Temptime_region = [];
TempFret12_region = [];
TempFret13_region = [];
TempFret23_region = [];
TempDonor_b_region = [];
TempDonor2_b_region = [];
TempAcceptor_b_region = [];
TempDonor_g_region = [];
TempDonor2_g_region = [];
TempAcceptor_g_region = [];
TempDonor_r_region = [];
TempDonor2_r_region = [];
TempAcceptor_r_region = [];

fileid_add_region = fopen(filename_add_region, 'r');
%fileid_add_region_data = fopen(filename_add_region_data, 'r');
fileid_add_region_data = -1;
if fileid_add_region~=-1 && fileid_add_region_data~=-1;
    fgetl(fileid_add_region);
    temp = fscanf(fileid_add_region, '%f %f %f  %f %f %f  %f %f %f  %f %f %f %f  %f %f  %f\n', [13 1]);
    
    dbackground_b=double(temp(1,1));
    d2background_b=double(temp(2,1));
    abackground_b=double(temp(3,1));
    
    dbackground_g=double(temp(4,1));
    d2background_g=double(temp(5,1));
    abackground_g=double(temp(6,1));
    
    dbackground_r=double(temp(7,1));
    d2background_r=double(temp(8,1));
    abackground_r=double(temp(9,1));
    
    leakage12=double(temp(10,1));   %0.11
    leakage21=double(temp(11,1));
    leakage13=double(temp(12,1));   %0.013
    leakage23=double(temp(13,1));   %0.12
    
    gamma12=double(temp(14,1));  %1
    gamma13=double(temp(15,1));  %2.36
    
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
time_length_each = floor(time_length/3);
DonorRawData_0 = zeros(NumberofPeaks, time_length_each, 'double');
Donor2RawData_0 = zeros(NumberofPeaks, time_length_each, 'double');
AcceptorRawData_0 = zeros(NumberofPeaks, time_length_each, 'double');
DonorRawData_1 = zeros(NumberofPeaks, time_length_each, 'double');
Donor2RawData_1 = zeros(NumberofPeaks, time_length_each, 'double');
AcceptorRawData_1 = zeros(NumberofPeaks, time_length_each, 'double');
DonorRawData_2 = zeros(NumberofPeaks, time_length_each, 'double');
Donor2RawData_2 = zeros(NumberofPeaks, time_length_each, 'double');
AcceptorRawData_2 = zeros(NumberofPeaks, time_length_each, 'double');

binlength = int32(time_length/binwidth-1);
bintime = zeros(binlength, 1, 'double');
binEraw = zeros(binlength, 1, 'double');
binEcorrect = zeros(binlength, 1, 'double');


for m=1:binlength
    bintime(m) = double(m-1)*(binwidth*timeunit);
end

% Rawdata�� ����� trace�����
if ColorNumber == 2
    for i=1:NumberofPeaks
        for j=1:time_length_each
            DonorRawData_1(i,j) = Data(i*2-1,j*2-1);
            Donor2RawData_1(i,j) = Data(i*2,j*2-1);
            %AcceptorRawData_1(i,j) = Data(i*2,j*2-1);
            AcceptorRawData_1(i,j) = 0;
            DonorRawData_2(i,j) = Data(i*2-1,j*2);
            Donor2RawData_2(i,j) = Data(i*2,j*2);
            %AcceptorRawData_2(i,j) = Data(i*2,j*2);
            AcceptorRawData_2(i,j) = 0;
        end
    end
end
if ColorNumber == 3
    for i=1:NumberofPeaks
        for j=1:time_length_each
            if LaserOrderChange == 'y' % Edit this part when the order of laser is wrong.
                DonorRawData_1(i,j) = Data(i*3-2,j*3-2);
                Donor2RawData_1(i,j) = Data(i*3-1,j*3-2);
                AcceptorRawData_1(i,j) = Data(i*3,j*3-2);
                DonorRawData_2(i,j) = Data(i*3-2,j*3-1);
                Donor2RawData_2(i,j) = Data(i*3-1,j*3-1);
                AcceptorRawData_2(i,j) = Data(i*3,j*3-1);
                DonorRawData_0(i,j) = Data(i*3-2,j*3);
                Donor2RawData_0(i,j) = Data(i*3-1,j*3);
                AcceptorRawData_0(i,j) = Data(i*3,j*3);
            else
                DonorRawData_0(i,j) = Data(i*3-2,j*3-2);
                Donor2RawData_0(i,j) = Data(i*3-1,j*3-2);
                AcceptorRawData_0(i,j) = Data(i*3,j*3-2);
                DonorRawData_1(i,j) = Data(i*3-2,j*3-1);
                Donor2RawData_1(i,j) = Data(i*3-1,j*3-1);
                AcceptorRawData_1(i,j) = Data(i*3,j*3-1);
                DonorRawData_2(i,j) = Data(i*3-2,j*3);
                Donor2RawData_2(i,j) = Data(i*3-1,j*3);
                AcceptorRawData_2(i,j) = Data(i*3,j*3);
            end
        end
    end
end

clear Data;



%% calculate, plot and save average traces %%
%MJ edited
% �� trace�� average�� ������, ù��° laser excitation�� Cy3 channel(����)�� �ι�° laser
%excitation�� cy3 channel�� ���� ��, ū ���� trace�� green laser excitation���� ���� ����
%trace�� red laser excitation���� ����
DonorRawData_b = DonorRawData_1;
Donor2RawData_b = Donor2RawData_1;
AcceptorRawData_b = AcceptorRawData_1;
DonorRawData_g = DonorRawData_2;
Donor2RawData_g = Donor2RawData_2;
AcceptorRawData_g = AcceptorRawData_2;
DonorRawData_r = DonorRawData_0;
Donor2RawData_r = Donor2RawData_0;
AcceptorRawData_r = AcceptorRawData_0;

% if sum(sum(DonorRawData_1, 1)) > sum(sum(DonorRawData_2, 1))
%     DonorRawData_g = DonorRawData_1;
%     Donor2RawData_g = Donor2RawData_1;
%     AcceptorRawData_g = AcceptorRawData_1;
%     DonorRawData_r = DonorRawData_2;
%     Donor2RawData_r = Donor2RawData_2;
%     AcceptorRawData_r = AcceptorRawData_2;
%     color_order = 0;
% else
%     DonorRawData_g = DonorRawData_2;
%     Donor2RawData_g = Donor2RawData_2;
%     AcceptorRawData_g = AcceptorRawData_2;
%     DonorRawData_r = DonorRawData_1;
%     Donor2RawData_r = Donor2RawData_1;
%     AcceptorRawData_r = AcceptorRawData_1;
% 	color_order = 1;
% end

DonorRawAverage_b = sum(DonorRawData_b, 1) / NumberofPeaks;
Donor2RawAverage_b = sum(Donor2RawData_b, 1) /NumberofPeaks;
AcceptorRawAverage_b = sum(AcceptorRawData_b, 1) / NumberofPeaks;
DonorRawAverage_g = sum(DonorRawData_g, 1) / NumberofPeaks;
Donor2RawAverage_g = sum(Donor2RawData_g, 1) /NumberofPeaks;
AcceptorRawAverage_g = sum(AcceptorRawData_g, 1) / NumberofPeaks;
DonorRawAverage_r = sum(DonorRawData_r, 1) / NumberofPeaks;
Donor2RawAverage_r = sum(Donor2RawData_r, 1) /NumberofPeaks;
AcceptorRawAverage_r = sum(AcceptorRawData_r, 1) / NumberofPeaks;


clear DonorRawData_0;
clear Donor2RawData_0;
clear AcceptorRawData_0;
clear DonorRawData_1;
clear Donor2RawData_1;
clear AcceptorRawData_1;
clear DonorRawData_2;
clear Donor2RawData_2;
clear AcceptorRawData_2;

figure('Name','Raw Data Ensemble Average');
hdl1 = gcf;

%Green excitation�� signal average
subplot(3,1,1);
plot(time_b, DonorRawAverage_b - dbackground_b, 'g', time_b, Donor2RawAverage_b - d2background_b, 'r', time_b, AcceptorRawAverage_b - abackground_b, 'b');
title('Raw data Average donor, donor2 and acceptor signal in Green excitation');
zoom on;

%Green excitation�� signal average
subplot(3,1,2);
plot(time_g, DonorRawAverage_g - dbackground_g, 'g', time_g, Donor2RawAverage_g - d2background_g, 'r', time_g, AcceptorRawAverage_g - abackground_g, 'b');
title('Raw data Average donor, donor2 and acceptor signal in Green excitation');
zoom on;

%Red excitation�� signal average
subplot(3,1,3);
plot(time_r, DonorRawAverage_r - dbackground_r, 'g', time_r, Donor2RawAverage_r - d2background_r, 'r', time_r, AcceptorRawAverage_r - abackground_r, 'b');
title('Raw data Average donor, donor2 and acceptor signal in Red excitation');
zoom on;

if Is_Avg_and_E_save =='n'
    AverageOutput_b = [time_b' DonorRawAverage_b' Donor2RawAverage_b' AcceptorRawAverage_b'];
    AverageOutput_g = [time_g' DonorRawAverage_g' Donor2RawAverage_g' AcceptorRawAverage_g'];
    AverageOutput_r = [time_r' DonorRawAverage_r' Donor2RawAverage_r' AcceptorRawAverage_r'];
    save([filename_head '_avg_b.dat'], 'AverageOutput_b', '-ascii');
    save([filename_head '_avg_g.dat'], 'AverageOutput_g', '-ascii');
    save([filename_head '_avg_r.dat'], 'AverageOutput_r', '-ascii');
end

%% calculate E level from the first 10 points and plot histograms of E level and total intensity. Also save the same info
tempDonor_b = reshape(DonorRawData_b(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempDonor2_b = reshape(Donor2RawData_b(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempAcceptor_b = reshape(AcceptorRawData_b(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempDonor_g = reshape(DonorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempDonor2_g = reshape(Donor2RawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempAcceptor_g = reshape(AcceptorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempDonor_r = reshape(DonorRawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempDonor2_r = reshape(Donor2RawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
tempAcceptor_r = reshape(AcceptorRawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);

if DyeType == 'cy235'
    EachTotal_23_g = tempDonor2_g + tempAcceptor_g;
    EachTotal_23_g = (EachTotal_23_g~=0).*EachTotal_23_g + (EachTotal_23_g==0)*1;	% remove zeros
    EachTotal_123_b = tempDonor_b + tempDonor2_b + tempAcceptor_b;
    EachTotal_123_b = (EachTotal_123_b~=0).*EachTotal_123_b + (EachTotal_123_b==0)*1;	% remove zeros
    
    E_level_23 = tempAcceptor_g./EachTotal_23_g;
    E_level_12 = tempDonor2_b./((1-E_level_23).*tempDonor_b + tempDonor2_b);
    E_level_13 = (tempAcceptor_b - E_level_23.*(tempDonor2_b + tempAcceptor_b))./(tempDonor_b + tempAcceptor_b - E_level_23.*(EachTotal_123_b));
    if Is_Avg_and_E_save =='y'
        E_level_output_b = [E_level_b EachTotal_b];
        E_level_output_g = [E_level_g EachTotal_g];
        save([filename_head '_elevel_10p_b.dat'],'E_level_output_b','-ascii');
        save([filename_head '_elevel_10p_g.dat'],'E_level_output_g','-ascii');
    end
end


figure('Name','Raw data analysis');
hdl2 = gcf;

subplot(3,4,1); % 2*3 figure_ upper left (the last number shows the location of figure)
hist(E_level_12,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw FRET_12 histogram' ]);
zoom on;

subplot(3,4,2); % 2*3 figure_ upper left (the last number shows the location of figure)
hist(E_level_13,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw FRET_13 histogram' ]);
zoom on;

subplot(3,4,3); % 2*3 figure_ upper left (the last number shows the location of figure)
hist(E_level_23,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw FRET_23 histogram' ]);
zoom on;

subplot(3,4,4);
hist(EachTotal_123_b,-100:50:4000);
temp=axis;
temp(1)=-100;
temp(2)=4000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Total intensity histogram' ]);
zoom on;

subplot(2,2,3);
plot(E_level_12, EachTotal_123_b,'b+', 'MarkerSize', 2);
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
temp(3)=-100;
temp(4)=4000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Total Intensity_123 vs. FRET' ]);
zoom on;

DonorFirstData_g = reshape(DonorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
Donor2FirstData_g = reshape(Donor2RawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
AcceptorFirstData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
size(DonorRawData_g)
time_length_each
DonorLastData_g = reshape(DonorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
Donor2LastData_g = reshape(Donor2RawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
AcceptorLastData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
DonorFirstData_r = reshape(DonorRawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
Donor2FirstData_r = reshape(Donor2RawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
AcceptorFirstData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
DonorLastData_r = reshape(DonorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
Donor2LastData_r = reshape(Donor2RawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
AcceptorLastData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);


subplot(5,4,11);
hist(DonorFirstData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Donor Intensity histogram']);
zoom on;

subplot(5,4,12);
hist(DonorLastData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Donor Intensity histogram']);
zoom on;

subplot(5,4,15);
hist(Donor2FirstData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Donor2 Intensity histogram']);
zoom on;

subplot(5,4,16);
hist(Donor2LastData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Donor2 Intensity histogram']);
zoom on;

subplot(5,4,19);
hist(AcceptorFirstData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Acceptor Intensity histogram']);
zoom on;

subplot(5,4,20);
hist(AcceptorLastData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Acceptor Intensity histogram']);
zoom on;


figure('Name','Raw data analysis additional part');
hdl_add = gcf;

subplot(2,2,1);
plot(E_level_12, EachTotal_123_b,'b+', 'MarkerSize', 2);
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
temp(3)=-100;
temp(4)=4000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Total Intensity_123 vs. FRET' ]);
zoom on;


subplot(2,2,3);
plot(Donor2FirstData_g, DonorFirstData_g,'b+', 'MarkerSize', 2);
temp=axis;
temp(1)=-100;
temp(2)=3000;
temp(3)=-100;
temp(4)=3000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Donor2 Intensity vs. Donor Intensity ' ]);
zoom on;


DonorFirstData_g_after = reshape(DonorRawData_g(1:NumberofPeaks, 2:(FirstNumber+1)), NumberofPeaks*FirstNumber, 1);
DonorFirstData_g_before = reshape(DonorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
Donor2FirstData_g_after = reshape(Donor2RawData_g(1:NumberofPeaks, 2:(FirstNumber+1)), NumberofPeaks*FirstNumber, 1);
Donor2FirstData_g_before = reshape(Donor2RawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);

E_level_12_after = Donor2FirstData_g_after./(Donor2FirstData_g_after + DonorFirstData_g_after);
E_level_12_before = Donor2FirstData_g_before./(Donor2FirstData_g_before + DonorFirstData_g_before);

E_level_12_onestep = [E_level_12(end)' E_level_12(1:(end-1))'];
subplot(2,2,2);
plot(E_level_12_before, E_level_12_after,'b+', 'MarkerSize', 2);
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
temp(3)=-0.1;
temp(4)=1.1;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Steps' ]);
zoom on;


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
DT1a=[];DT2a=[];DT3a=[];
DT1d=[];DT2d=[];DT3d=[];
DT1f=[];DT2f=[];DT3f=[];

scrsz = get(0,'ScreenSize');
figure('Name','Trace analysis','OuterPosition',[scrsz(4)/2 0.05*scrsz(4) scrsz(4)+300 0.95*scrsz(4)]);
hdl_trace=gcf;

i=0;
history_n = 0;
history = zeros(1000, 1, 'int16');

% Display traces

while i < NumberofPeaks
    i = int32(int32(i) + 1);
    
    again=1;
    while ((again==1) && (i <= NumberofPeaks))
        if ColorNumber == 3
            if DyeType == 'cy235'
                DonorCorrect_b = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_b(i,:) - dbackground_b) - leakage21 * (Donor2RawData_g(i,:) - d2background_b));
                Donor2Correct_b = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_b(i,:) - d2background_b) - leakage12 * (DonorRawData_b(i,:) - dbackground_b));
                AcceptorCorrect_b = gamma13 * ((AcceptorRawData_b(i,:) - abackground_b) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_b(i,:) - dbackground_b) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_b(i,:) -d2background_b));
                EachTotalCorrect_b = DonorCorrect_b + Donor2Correct_b + AcceptorCorrect_b;
                EachTotalCorrect_b = (EachTotalCorrect_b~=0).*EachTotalCorrect_b + (EachTotalCorrect_b==0)*1;	% remove zeros
                
                DonorCorrect_g = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_g(i,:) - dbackground_g) - leakage21 * (Donor2RawData_g(i,:) - d2background_g));
                Donor2Correct_g = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_g(i,:) - d2background_g) - leakage12 * (DonorRawData_g(i,:) - dbackground_g));
                AcceptorCorrect_g = gamma13 * ((AcceptorRawData_g(i,:) - abackground_g) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_g(i,:) - dbackground_g) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_g(i,:) -d2background_g));
                AcceptorCorrect_g = AcceptorCorrect_g - direct * (Donor2Correct_g + AcceptorCorrect_g);
                EachTotalCorrect_g = Donor2Correct_g + AcceptorCorrect_g;
                EachTotalCorrect_g = (EachTotalCorrect_g~=0).*EachTotalCorrect_g + (EachTotalCorrect_g==0)*1;	% remove zeros
                
                DonorCorrect_r = DonorRawData_r(i,:);
                Donor2Correct_r = Donor2RawData_r(i,:);
                AcceptorCorrect_r = AcceptorRawData_r(i,:);
                EachTotalCorrect_r = AcceptorCorrect_r;
                EachTotalCorrect_r = (EachTotalCorrect_r~=0).*EachTotalCorrect_r + (EachTotalCorrect_r==0)*1;	% remove zeros
                
                Fret23 = AcceptorCorrect_g./EachTotalCorrect_g;
                Fret12 = Donor2Correct_b./((1-Fret23).*DonorCorrect_b + Donor2Correct_b);
                Fret13 = (AcceptorCorrect_b - Fret23.*(Donor2Correct_b + AcceptorCorrect_b))./(DonorCorrect_b + AcceptorCorrect_b - Fret23 .* (EachTotalCorrect_b));
            end
        end
        if ColorNumber == 2
            DonorCorrect_g = (DonorRawData_g(i,:) - dbackground_g) + leakage12 * (DonorRawData_g(i,:) - dbackground_g);
            Donor2Correct_g = gamma12 * ((Donor2RawData_g(i,:) - d2background_g) - leakage12 * (DonorRawData_g(i,:) - dbackground_g));
            AcceptorCorrect_g = 0 * DonorRawData_g(i,:); %xxx
            EachTotalCorrect_g = DonorCorrect_g + Donor2Correct_g + AcceptorCorrect_g;
            EachTotalCorrect_g = (EachTotalCorrect_g~=0).*EachTotalCorrect_g + (EachTotalCorrect_g==0)*1;	% remove zeros
            
            DonorCorrect_r = (DonorRawData_r(i,:) - dbackground_r) + leakage12 * (DonorRawData_r(i,:) - dbackground_r);
            Donor2Correct_r = gamma12 * ((Donor2RawData_r(i,:) - d2background_r) - leakage12 * (DonorRawData_r(i,:) - dbackground_r));
            AcceptorCorrect_r = 0 * DonorCorrect_r;
            EachTotalCorrect_r = DonorCorrect_r + Donor2Correct_r;
            EachTotalCorrect_r = (EachTotalCorrect_r~=0).*EachTotalCorrect_r + (EachTotalCorrect_r==0)*1;	% remove zeros
            
            Fret23 = 0 * DonorCorrect_g;
            Fret12 = Donor2Correct_g./EachTotalCorrect_g;
            Fret13 = 0 * DonorCorrect_g;
        end
        
        %for j=2:time_length_each
        %    if (Fret13(j) < -0.3 | Fret13(j) > 1.1 | Fret12(j) < -0.3 | Fret23(j) < -0.3 | Fret12(j) > 1.1 | Fret23(j) > 1.1)
        %        Fret12(j) = Fret12(j-1);
        %        Fret13(j) = Fret13(j-1);
        %        Fret23(j) = Fret23(j-1);
        %    end
        %end
        
        for j=2:time_length_each
            if (Fret13(j) < -0.3)
                Fret13(j) = 0;
            elseif (Fret13(j) > 1.1)
                Fret13(j) = 1;
            end
            if(Fret12(j) < -0.3)
                Fret12(j) = 0;
            elseif(Fret12(j) > 1.1)
                Fret12(j) = 1;
            end
            if(Fret23(j) < -0.3)
                Fret23(j) = 0;
            elseif(Fret23(j) > 1.1)
                Fret23(j) = 1;
            end
        end
        
        
        %FretEcCorrect_g = AcceptorCorrect_g./EachTotalCorrect_g;
        
        %EachTotalCorrect_r = DonorCorrect_r + Donor2Correct_r + AcceptorCorrect_r;
        %EachTotalCorrect_r = (EachTotalCorrect_r~=0).*EachTotalCorrect_r + (EachTotalCorrect_r==0)*1;	% remove zeros
        %FretEcCorrect_r = AcceptorCorrect_r./EachTotalCorrect_r;
        %		if DoseBinningNeed == 'y'
        %			for m=1:binlength
        %				binEcorrect_g(m) = sum(FretEcCorrect_g((double(m-1)*binwidth+1):(binwidth*m) ) )/binwidth;
        %               binEcorrect_r(m) = sum(FretEcCorrect_r((double(m-1)*binwidth+1):(binwidth*m) ) )/binwidth;
        %			end
        %		end
        
        
        %FretEc_filter = 0.3;
        %MeanRange = 100;
        %Acceptor_filter = 100;
        %Total_filter = 500;
        %if DoesFilterNeed ~= 'y' | ( mean(FretEcCorrect_g(1:MeanRange)) > FretEc_filter && mean(AcceptorCorrect_g(1:MeanRange)) > Acceptor_filter && mean(EachTotalCorrect_g(1:MeanRange)) > Total_filter )
        %	i = i - 1;
        again=0;
        %end
        %i = i + 1 ;
    end
    
    % Trace window
    % blue laser excitation corrected trace
    figure(hdl_trace);
    subplot('position',[0.1 0.84 0.8 0.10]); %first trace : Blue excitation 3-color corrected trace
    plot(time_b, DonorCorrect_b, 'g', time_b, Donor2Correct_b, 'r', time_b, AcceptorCorrect_b, 'm');
    hold on
    plot(time_b, DonorCorrect_b + Donor2Correct_b + AcceptorCorrect_b + 300, 'k');
    hold off
    temp=axis;
    temp(3)=BottomLimit_b;
    temp(4)=UpperLimit_b;
    grid on;
    axis(temp);
    %title(['                    Molecule ' num2str(i) '  / ' num2str(NumberofPeaks) '              file: ' filename_head ]); %'   cor.: ' num2str(dbackground) '  ' num2str(abackground) '  ' num2str(d2background) ' ' num2str(leakage12) '  ' num2str(leakage13) ' ' num2str(leakage23) ' ' num2str(gamma12) ' ' num2str(gamma13) '   filter: ' DoesFilterNeed ' Fret:' num2str(FretEc_filter) ' Mean:' num2str(MeanRange) ' Acceptor:' num2str(Acceptor_filter) ' Total:' num2str(Total_filter) ]);
    title(['ALEX Green Laser Molecule' num2str(i) '  / ' num2str(NumberofPeaks) '       file:' filename_head  ' cor.: ' num2str(dbackground_b) '  ' num2str(d2background_b) '  ' num2str(abackground_b) '  ' num2str(dbackground_g) '  ' num2str(d2background_g) '  ' num2str(abackground_g) '  ' num2str(dbackground_r) '  ' num2str(d2background_r) '  ' num2str(abackground_r) '  ' num2str(leakage12) '  ' num2str(leakage21) '  ' num2str(leakage13) '  ' num2str(leakage23) '  ' num2str(gamma12) '  ' num2str(gamma13) ]);
    zoom on;
    
    region_indice = find(Temp_i == i);
    number_region = size(region_indice);
    if number_region ~= 0
        FirstSelectX=[(Tempfirstpoint(region_indice)*timeunit)' (Tempfirstpoint(region_indice)*timeunit)'];
        FirstSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
        LastSelectX=[(Templastpoint(region_indice)*timeunit)' (Templastpoint(region_indice)*timeunit)'];
        LastSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
    end
    
    % green laser excitation corrected trace
    subplot('position',[0.1 0.72 0.8 0.10]); %first trace : Green excitation 3-color corrected trace
    plot(time_g, DonorCorrect_g, 'g', time_g, Donor2Correct_g, 'r', time_g, AcceptorCorrect_g, 'm');
    temp=axis;
    temp(3)=BottomLimit_g;
    temp(4)=UpperLimit_g;
    grid on;
    axis(temp);
    %title(['                    Molecule ' num2str(i) '  / ' num2str(NumberofPeaks) '              file: ' filename_head ]); %'   cor.: ' num2str(dbackground) '  ' num2str(abackground) '  ' num2str(d2background) ' ' num2str(leakage12) '  ' num2str(leakage13) ' ' num2str(leakage23) ' ' num2str(gamma12) ' ' num2str(gamma13) '   filter: ' DoesFilterNeed ' Fret:' num2str(FretEc_filter) ' Mean:' num2str(MeanRange) ' Acceptor:' num2str(Acceptor_filter) ' Total:' num2str(Total_filter) ]);
    title(['ALEX Red Laser Molecule' num2str(i) ' / ' num2str(NumberofPeaks) ' file:' filename_head  ' cor.: ' num2str(dbackground_g) '  ' num2str(d2background_g) '  ' num2str(abackground_g) '  ' num2str(dbackground_r) '  ' num2str(d2background_r) '  ' num2str(abackground_r) '  ' num2str(leakage12) '  ' num2str(leakage21) '  ' num2str(leakage13) '  ' num2str(leakage23) '  ' num2str(gamma12) '  ' num2str(gamma13) ]);
    zoom on;
    
    region_indice = find(Temp_i == i);
    number_region = size(region_indice);
    if number_region ~= 0
        FirstSelectX=[(Tempfirstpoint(region_indice)*timeunit)' (Tempfirstpoint(region_indice)*timeunit)'];
        FirstSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
        LastSelectX=[(Templastpoint(region_indice)*timeunit)' (Templastpoint(region_indice)*timeunit)'];
        LastSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
    end
    
    % red laser excitation corrected trace
    subplot('position',[0.1 0.60 0.8 0.10]);
    plot(time_g, DonorCorrect_r, 'g', time_g, Donor2Correct_r, 'r', time_g, AcceptorCorrect_r, 'm');
    temp=axis;
    temp(3)=BottomLimit_r;
    temp(4)=UpperLimit_r;
    grid on;
    axis(temp);
    title(['ALEX 750 Laser Molecule' num2str(i) ' / ' num2str(NumberofPeaks) ' timeunit: ' num2str(timeunit) '  gain: ' num2str(gain) '  scaler: ' num2str(scaler)  '  spot diameter: ' num2str(SpotDiameter) '    ' date]);
    zoom on;
    
    
    subplot('position',[0.93 0.42 0.03 0.15]);
    x = -0.1:0.02:1.1;
    [hX,hN]=hist(Fret12,x);
    barh(hN,hX,'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    axis on;
    zoom on;
    title('FRET12')
    
    region_indice = find(Temp_i == i);
    number_region = sum(size(region_indice)) - 1;
    FirstSelectX=[];
    FirstSelectY=[];
    LastSelectX=[];
    LastSelectY=[];
    for j=1:number_region
        tempx=[Tempfirstpoint(region_indice(j))*2*timeunit Tempfirstpoint(region_indice(j))*2*timeunit];
        tempy=[ (2 - 4 * mod(j,2)) (4 * mod(j,2) - 2)];
        FirstSelectX=[FirstSelectX tempx];
        FirstSelectY=[FirstSelectY tempy];
        tempx=[Templastpoint(region_indice(j))*2*timeunit Templastpoint(region_indice(j))*2*timeunit];
        LastSelectX=[LastSelectX tempx];
        LastSelectY=[LastSelectY tempy];
    end
    
    subplot('position',[0.1 0.42 0.8 0.15]);
    %	FretEc=(1./(1+gamma*(donorcorrect(i,:)./acceptorcorrect(i,:))));
    hFretLine = plot(time_g, Fret12, FirstSelectX, FirstSelectY, LastSelectX, LastSelectY, bintime, binEcorrect, 'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    zoom on;
    title('FRET12')
    
    subplot('position',[0.1 0.22 0.8 0.15]);
    plot(time_g, Fret13, bintime, binEraw, 'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    zoom on;
    title('FRET13')
    
    subplot('position',[0.93 0.22 0.03 0.15]);
    x = -0.1:0.02:1.1;
    [hX,hN]=hist(Fret13,x);
    barh(hN,hX,'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    axis on;
    zoom on;
    title('FRET13')
    
    subplot('position', [0.1 0.02 0.8 0.15]);
    plot(time_g, Fret23, bintime, binEraw, 'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp)
    grid on;
    zoom on;
    title('FRET23')
    
    subplot('position',[0.93 0.02 0.03 0.15]);
    x = -0.1:0.02:1.1;
    [hX,hN]=hist(Fret23,x);
    barh(hN,hX,'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    axis on;
    zoom on;
    title('FRET23')
    
    if DoseMovieNeed == 'y'
        fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
        %startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + color_order*peak_height*peaks_total_width;
        Xpoint = 6;
        startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 0*peak_height*peaks_total_width;
        for j=1:peak_height
            fseek(fileid_movie, startpoint , 'bof');
            peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
            peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
            startpoint = startpoint + peaks_total_width;
        end
        fclose(fileid_movie);
        
        subplot('position',[0.93 0.80 0.05 0.18]);
        colormap(hot);
        image(peak);
        axis off;
        zoom on;
        title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
        
        fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
        %startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + (1-color_order)*peak_height*peaks_total_width;
        startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 1*peak_height*peaks_total_width;
        for j=1:peak_height
            fseek(fileid_movie, startpoint , 'bof');
            peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
            peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
            startpoint = startpoint + peaks_total_width;
        end
        fclose(fileid_movie);
        
        subplot('position',[0.93 0.60 0.05 0.18]);
        colormap(hot);
        image(peak);
        axis off;
        zoom on;
        title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
        
        fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
        %startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + (1-color_order)*peak_height*peaks_total_width;
        startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 2*peak_height*peaks_total_width;
        for j=1:peak_height
            fseek(fileid_movie, startpoint , 'bof');
            peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
            peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
            startpoint = startpoint + peaks_total_width;
        end
        fclose(fileid_movie);
        
        subplot('position',[0.93 0.40 0.05 0.18]);
        colormap(hot);
        image(peak);
        axis off;
        zoom on;
        title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
    end
    
    
    again=1;
    while again==1
        again=0;
        disp([num2str(i) ' (l=save select, u=save region, d=delete region, s=save region, h=histogram select, p=display movie)']);
        keyanswer =input('(t=terminate program, b=back, g=go, r=calculate gamma, o=subtract background, k=calculate leakage, i=calculate direction cy7 excitation, q=collect photobleaching time) : ','s');
        answer = sscanf(keyanswer, '%s %*s');
        numberofanswer = sscanf(keyanswer, '%*s %f');
        
        if answer=='o'
            
            [raw_x, ~] = ginput(2);
            x = round(raw_x / timeunit);
            dbackground_0_temp = mean(DonorRawData_0, x(1):x(2));
            d2background_0_temp = mean(Donor2RawData_0, x(1):x(2));
            abackground_0_temp = mean(AcceptorRawData_0, x(1):x(2));
            
            i = i - 1;
            
            
%             again=1;
%             [raw_x, raw_y] = ginput(2);
%             x = round(raw_x / (3*timeunit));
%             d1d1_bg = mean(DonorCorrect_b(x(1):x(2)));
%             d1d2_bg = mean(Donor2Correct_b(x(1):x(2)));
%             d1ac_bg = mean(AcceptorCorrect_b(x(1):x(2)));
%             
%             d2d1_bg = mean(DonorCorrect_g(x(1):x(2)));
%             d2d2_bg = mean(Donor2Correct_g(x(1):x(2)));
%             d2ac_bg = mean(AcceptorCorrect_g(x(1):x(2)));
%             
%             acd1_bg = mean(DonorCorrect_r(x(1):x(2)));
%             acd2_bg = mean(Donor2Correct_r(x(1):x(2)));
%             acac_bg = mean(AcceptorCorrect_r(x(1):x(2)));
%             
%             DonorCorrect_b = DonorCorrect_b - d1d1_bg;
%             Donor2Correct_b = Donor2Correct_b - d1d2_bg;
%             AcceptorCorrect_b = AcceptorCorrect_b - d1ac_bg;
%             
%             DonorCorrect_g = DonorCorrect_g - d2d1_bg;
%             Donor2Correct_g = Donor2Correct_g - d2d2_bg;
%             AcceptorCorrect_g = AcceptorCorrect_g - d2ac_bg;
%             
%             DonorCorrect_r = DonorCorrect_r - acd1_bg;
%             Donor2Correct_r = Donor2Correct_r - acd2_bg;
%             AcceptorCorrect_r = AcceptorCorrect_r - acac_bg;
%             
%             figure(hdl_trace);
%             subplot('position',[0.1 0.84 0.8 0.10]); %first trace : Blue excitation 3-color corrected trace
%             plot(time_b, DonorCorrect_b, 'g', time_b, Donor2Correct_b, 'r', time_b, AcceptorCorrect_b, 'm');
%             hold on
%             plot(time_b, DonorCorrect_b + Donor2Correct_b + AcceptorCorrect_b + 100, 'k');
%             hold off
%             temp=axis;
%             temp(3)=BottomLimit_b;
%             temp(4)=UpperLimit_b;
%             grid on;
%             axis(temp);
%             title(['ALEX Green Laser Molecule' num2str(i) '  / ' num2str(NumberofPeaks) '       file:' filename_head  ' cor.: ' num2str(dbackground_b) '  ' num2str(d2background_b) '  ' num2str(abackground_b) '  ' num2str(dbackground_g) '  ' num2str(d2background_g) '  ' num2str(abackground_g) '  ' num2str(dbackground_r) '  ' num2str(d2background_r) '  ' num2str(abackground_r) '  ' num2str(leakage12) '  ' num2str(leakage21) '  ' num2str(leakage13) '  ' num2str(leakage23) '  ' num2str(gamma12) '  ' num2str(gamma13) ]);
%             zoom on;
%             
%             region_indice = find(Temp_i == i);
%             number_region = size(region_indice);
%             if number_region ~= 0
%                 FirstSelectX=[(Tempfirstpoint(region_indice)*timeunit)' (Tempfirstpoint(region_indice)*timeunit)'];
%                 FirstSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
%                 LastSelectX=[(Templastpoint(region_indice)*timeunit)' (Templastpoint(region_indice)*timeunit)'];
%                 LastSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
%             end
%             
%             % green laser excitation corrected trace
%             subplot('position',[0.1 0.72 0.8 0.10]); %first trace : Green excitation 3-color corrected trace
%             plot(time_g, DonorCorrect_g, 'g', time_g, Donor2Correct_g, 'r', time_g, AcceptorCorrect_g, 'm');
%             temp=axis;
%             temp(3)=BottomLimit_g;
%             temp(4)=UpperLimit_g;
%             grid on;
%             axis(temp);
%             title(['ALEX Red Laser Molecule' num2str(i) ' / ' num2str(NumberofPeaks) ' file:' filename_head  ' cor.: ' num2str(dbackground_g) '  ' num2str(d2background_g) '  ' num2str(abackground_g) '  ' num2str(dbackground_r) '  ' num2str(d2background_r) '  ' num2str(abackground_r) '  ' num2str(leakage12) '  ' num2str(leakage21) '  ' num2str(leakage13) '  ' num2str(leakage23) '  ' num2str(gamma12) '  ' num2str(gamma13) ]);
%             zoom on;
%             
%             region_indice = find(Temp_i == i);
%             number_region = size(region_indice);
%             if number_region ~= 0
%                 FirstSelectX=[(Tempfirstpoint(region_indice)*timeunit)' (Tempfirstpoint(region_indice)*timeunit)'];
%                 FirstSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
%                 LastSelectX=[(Templastpoint(region_indice)*timeunit)' (Templastpoint(region_indice)*timeunit)'];
%                 LastSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
%             end
%             
%             % red laser excitation corrected trace
%             subplot('position',[0.1 0.60 0.8 0.10]);
%             plot(time_g, DonorCorrect_r, 'g', time_g, Donor2Correct_r, 'r', time_g, AcceptorCorrect_r, 'm');
%             temp=axis;
%             temp(3)=BottomLimit_r;
%             temp(4)=UpperLimit_r;
%             grid on;
%             axis(temp);
%             title(['ALEX 750 Laser Molecule' num2str(i) ' / ' num2str(NumberofPeaks) ' timeunit: ' num2str(timeunit) '  gain: ' num2str(gain) '  scaler: ' num2str(scaler)  '  spot diameter: ' num2str(SpotDiameter) '    ' date]);
%             zoom on;
%             
%             
%             subplot('position',[0.93 0.42 0.03 0.15]);
%             x = -0.1:0.02:1.1;
%             [hX,hN]=hist(Fret12,x);
%             barh(hN,hX,'k');
%             temp=axis;
%             temp(3)=-0.1;
%             temp(4)=1.1;
%             axis(temp);
%             grid on;
%             axis on;
%             zoom on;
%             title('FRET12')
%             
%             region_indice = find(Temp_i == i);
%             number_region = sum(size(region_indice)) - 1;
%             FirstSelectX=[];
%             FirstSelectY=[];
%             LastSelectX=[];
%             LastSelectY=[];
%             for j=1:number_region
%                 tempx=[Tempfirstpoint(region_indice(j))*2*timeunit Tempfirstpoint(region_indice(j))*2*timeunit];
%                 tempy=[ (2 - 4 * mod(j,2)) (4 * mod(j,2) - 2)];
%                 FirstSelectX=[FirstSelectX tempx];
%                 FirstSelectY=[FirstSelectY tempy];
%                 tempx=[Templastpoint(region_indice(j))*2*timeunit Templastpoint(region_indice(j))*2*timeunit];
%                 LastSelectX=[LastSelectX tempx];
%                 LastSelectY=[LastSelectY tempy];
%             end
%             
%             subplot('position',[0.1 0.42 0.8 0.15]);
%             %	FretEc=(1./(1+gamma*(donorcorrect(i,:)./acceptorcorrect(i,:))));
%             hFretLine = plot(time_g, Fret12, FirstSelectX, FirstSelectY, LastSelectX, LastSelectY, bintime, binEcorrect, 'k');
%             temp=axis;
%             temp(3)=-0.1;
%             temp(4)=1.1;
%             axis(temp);
%             grid on;
%             zoom on;
%             title('FRET12')
%             
%             subplot('position',[0.1 0.22 0.8 0.15]);
%             plot(time_g, Fret13, bintime, binEraw, 'k');
%             temp=axis;
%             temp(3)=-0.1;
%             temp(4)=1.1;
%             axis(temp);
%             grid on;
%             zoom on;
%             title('FRET13')
%             
%             subplot('position',[0.93 0.22 0.03 0.15]);
%             x = -0.1:0.02:1.1;
%             [hX,hN]=hist(Fret13,x);
%             barh(hN,hX,'k');
%             temp=axis;
%             temp(3)=-0.1;
%             temp(4)=1.1;
%             axis(temp);
%             grid on;
%             axis on;
%             zoom on;
%             title('FRET13')
%             
%             subplot('position', [0.1 0.02 0.8 0.15]);
%             plot(time_g, Fret23, bintime, binEraw, 'k');
%             temp=axis;
%             temp(3)=-0.1;
%             temp(4)=1.1;
%             axis(temp)
%             grid on;
%             zoom on;
%             title('FRET23')
%             
%             subplot('position',[0.93 0.02 0.03 0.15]);
%             x = -0.1:0.02:1.1;
%             [hX,hN]=hist(Fret23,x);
%             barh(hN,hX,'k');
%             temp=axis;
%             temp(3)=-0.1;
%             temp(4)=1.1;
%             axis(temp);
%             grid on;
%             axis on;
%             zoom on;
%             title('FRET23')
            
            if DoseMovieNeed == 'y'
                fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
                %startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + color_order*peak_height*peaks_total_width;
                Xpoint = 6;
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 0*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
                end
                fclose(fileid_movie);
                
                subplot('position',[0.93 0.80 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                zoom on;
                title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
                
                fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
                %startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + (1-color_order)*peak_height*peaks_total_width;
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 1*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
                end
                fclose(fileid_movie);
                
                subplot('position',[0.93 0.60 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                zoom on;
                title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
                
                fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
                %startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + (1-color_order)*peak_height*peaks_total_width;
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 2*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
                end
                fclose(fileid_movie);
                
                subplot('position',[0.93 0.40 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                zoom on;
                title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
            end
            
        end
        
        if answer=='p'
            again=1;
            [Xc,Yc] = ginput(1);
            Xpoint = round(Xc(1)/timeunit);
            
            subplot('position',[0.1 0.42 0.8 0.15]);
            % fretEc=(1./(1+gamma*(donorcorrect(i,:)./acceptorcorrect(i,:))));
            SelectX=[Xc Xc];
            SelectY=[-2 +2];
            plot(time_g, Fret12, SelectX, SelectY, 'k');
            temp=axis;
            temp(3)=-0.1;
            temp(4)=1.1;
            axis(temp);
            grid on;
            zoom on;
            
            if DoseMovieNeed == 'y'
                fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
                
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 0*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint, 'bof');
                    peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
                end
                fclose(fileid_movie);
                
                subplot('position',[0.93 0.80 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                zoom on;
                title(['frame: ' num2str(3*(uint32(Xpoint/3)-1)) '  time ' num2str(Xpoint*timeunit)]);
                
                fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 1*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
                end
                fclose(fileid_movie);
                
                subplot('position',[0.93 0.60 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                zoom on;
                title(['frame: ' num2str(3*(uint32(Xpoint/3)-1)+1) '  time ' num2str(Xpoint*timeunit)]);
                
                fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_3alex_movies'], 'r');
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 3*(uint32(Xpoint/3)-1)*peak_height*peaks_total_width + 2*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
                end
                fclose(fileid_movie);
                
                subplot('position',[0.93 0.40 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                zoom on;
                title(['frame: ' num2str(3*(uint32(Xpoint/3)-1)+2) '  time ' num2str(Xpoint*timeunit)]);
            end
        end
        
        if answer=='h'
            again=1;
            [Xc,Yc] = ginput(1);
            firstpoint = round(Xc(1)/(2*timeunit));
            
            subplot('position',[0.1 0.42 0.8 0.15]);
            %fretEc=(1./(1+gamma*(donorcorrect(i,:)./acceptorcorrect(i,:))));
            
            FirstSelectX=[Xc Xc];
            FirstSelectY=[-2 +2];
            plot(time_g, Fret12, FirstSelectX, FirstSelectY, bintime, binEcorrect, 'k');
            temp=axis;
            temp(3)=-0.1;
            temp(4)=1.1;
            axis(temp);
            grid on;
            zoom on;
            
            [Xc,Yc] = ginput(1);
            lastpoint = round(Xc(1)/(2*timeunit));
            
            subplot('position',[0.1 0.42 0.8 0.15]);
            LastSelectX=[Xc Xc];
            LastSelectY=[-2 +2];
            plot(time_g, Fret12, FirstSelectX, FirstSelectY, LastSelectX, LastSelectY, bintime, binEcorrect, 'k');
            temp=axis;
            temp(3)=-0.1;
            temp(4)=1.1;
            axis(temp);
            grid on;
            zoom on;
            
            if firstpoint>lastpoint
                temp=lastpoint;
                lastpoint=firstpoint;
                firstpoint=temp;
            end
            
            disp(firstpoint);
            disp(lastpoint);
            subplot('position',[0.93 0.42 0.06 0.15]);
            x = -0.1:0.02:1.1;
            [hX,hN]=hist(Fret12(firstpoint:lastpoint),x);
            barh(hN,hX,'k');
            temp=axis;
            temp(3)=-0.1;
            temp(4)=1.1;
            axis(temp);
            grid on;
            axis on;
            zoom on;
        end
        
        if answer=='u' | answer=='d'
            filename_add_region=[ filename_head '_region.dat' ];
            filename_add_region_data=[ filename_head '_region_data.dat' ];
            if answer=='u'
                [Xc,Yc,buttonc] = ginput(2); %ginput means graphicalinput from a mouse or cursor
                firstpoint = round(Xc(1)/(2*timeunit)); % fp means the Akaike final prediction for estimate
                lastpoint = round(Xc(2)/(2*timeunit)); % get database column privilige
                
                Temp_i = [Temp_i i];
                Tempfirstpoint = [Tempfirstpoint firstpoint];
                Templastpoint = [Templastpoint lastpoint];
                Templength = [Templength (lastpoint-firstpoint+1)];
            end
            
            if answer=='d'
                del_indice = find(Temp_i == i);
                sizeofarray = sum(size(del_indice)) - 1;
                if sizeofarray ~= 0
                    disp(Tempfirstpoint(del_indice));
                    temp=str2num(input(['number to delete (' int2str(1:sizeofarray) ')  :'],'s'));
                    del_indice = del_indice(temp);
                    sizeofarray = sum(size(Temp_i)) - 1;
                    select_indice = [ 1:(del_indice-1) (del_indice+1):sizeofarray ];
                    Temp_i = Temp_i( select_indice );
                    Tempfirstpoint = Tempfirstpoint( select_indice );
                    Templastpoint = Templastpoint( select_indice );
                    Templength = Templength( select_indice );
                end
            end
            
            Temptime_region = [];
            TempFret12_region = [];
            TempFret13_region = [];
            TempFret23_region = [];
            TempDonor_g_region = [];
            TempDonor2_g_region = [];
            TempDonor_r_region = [];
            TempDonor2_r_region = [];
            TempAcceptor_g_region = [];
            TempAcceptor_r_region = [];
            
            sizeofarray = sum(size(Temp_i)) - 1;
            size(Tempfirstpoint)
            for j=1:sizeofarray
                if ColorNumber == 2
                    temp_DonorCorrect_g = (DonorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_g) + leakage12 * (DonorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_g);
                    temp_Donor2Correct_g = gamma12 * ((Donor2RawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - d2background_g) - leakage12 * (DonorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_g));
                    temp_AcceptorCorrect_g = 0 * DonorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j));
                    temp_EachTotalCorrect_g = temp_DonorCorrect_g + temp_Donor2Correct_g;
                    temp_EachTotalCorrect_g = (temp_EachTotalCorrect_g~=0).*temp_EachTotalCorrect_g + (temp_EachTotalCorrect_g==0)*1;	% remove zeros
                    
                    temp_DonorCorrect_r = (DonorRawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_r) + leakage12 * (DonorRawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_r);
                    temp_Donor2Correct_r = gamma12 * ((Donor2RawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - d2background_r) - leakage12 * (DonorRawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_r));
                    temp_AcceptorCorrect_r = 0 * temp_DonorCorrect_r;
                    temp_EachTotalCorrect_r = temp_DonorCorrect_r + temp_Donor2Correct_r;
                    temp_EachTotalCorrect_r = (temp_EachTotalCorrect_r~=0).*temp_EachTotalCorrect_r + (temp_EachTotalCorrect_r==0)*1;	% remove zeros
                    
                    temp_EachTotalCorrect_g = temp_DonorCorrect_g + temp_Donor2Correct_g + temp_AcceptorCorrect_g;
                    temp_EachTotalCorrect_g = (temp_EachTotalCorrect_g~=0).*temp_EachTotalCorrect_g + (temp_EachTotalCorrect_g==0)*1;	% remove zeros
                    temp_Fret23 = 0 * temp_DonorCorrect_g;
                    temp_Fret12 = temp_Donor2Correct_g./temp_EachTotalCorrect_g;
                    temp_Fret13 = 0 * temp_DonorCorrect_g;
                end
                if ColorNumber == 3
                    temp_DonorCorrect_g = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_g) - leakage21 * (Donor2RawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - d2background_g));
                    temp_Donor2Correct_g = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - d2background_g) - leakage12 * (DonorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_g));
                    temp_AcceptorCorrect_g = gamma13 * ((AcceptorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - abackground_g) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_g) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_g(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) -d2background_g));
                    temp_EachTotalCorrect_g = temp_DonorCorrect_g + temp_Donor2Correct_g + temp_AcceptorCorrect_g;
                    temp_EachTotalCorrect_g = (temp_EachTotalCorrect_g~=0).*temp_EachTotalCorrect_g + (temp_EachTotalCorrect_g==0)*1;	% remove zeros
                    
                    temp_DonorCorrect_r = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_r) - leakage21 * (Donor2RawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - d2background_r));
                    temp_Donor2Correct_r = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - d2background_r) - leakage12 * (DonorRawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_r));
                    temp_AcceptorCorrect_r = gamma13 * ((AcceptorRawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - abackground_r) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) - dbackground_r) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_r(Temp_i(j),Tempfirstpoint(j):Templastpoint(j)) -d2background_r));
                    temp_AcceptorCorrect_r = temp_AcceptorCorrect_r - direct * (temp_Donor2Correct_r + temp_AcceptorCorrect_r);
                    temp_EachTotalCorrect_r = temp_Donor2Correct_r + temp_AcceptorCorrect_r;
                    temp_EachTotalCorrect_r = (temp_EachTotalCorrect_r~=0).*temp_EachTotalCorrect_r + (temp_EachTotalCorrect_r==0)*1;	% remove zeros
                    
                    temp_Fret23 = temp_AcceptorCorrect_r./temp_EachTotalCorrect_r;
                    temp_Fret12 = temp_Donor2Correct_g./((1-temp_Fret23).*temp_DonorCorrect_g + temp_Donor2Correct_g);
                    temp_Fret13 = (temp_AcceptorCorrect_g - temp_Fret23.*(temp_Donor2Correct_g + temp_AcceptorCorrect_g))./(temp_DonorCorrect_g + temp_AcceptorCorrect_g - temp_Fret23 .* (temp_EachTotalCorrect_g));
                end
                
                Temptime_region = [Temptime_region time(Tempfirstpoint(j):Templastpoint(j)) ];
                TempFret12_region = [TempFret12_region temp_Fret12 ];
                TempFret13_region = [TempFret13_region temp_Fret13 ];
                TempFret23_region = [TempFret23_region temp_Fret23 ];
                TempDonor_g_region = [TempDonor_g_region temp_DonorCorrect_g ];
                TempDonor2_g_region = [TempDonor2_g_region temp_Donor2Correct_g ];
                TempDonor_r_region = [TempDonor_r_region temp_DonorCorrect_r ];
                TempDonor2_r_region = [TempDonor2_r_region temp_Donor2Correct_r ];
                TempAcceptor_g_region = [TempAcceptor_g_region temp_AcceptorCorrect_g ];
                TempAcceptor_r_region = [TempAcceptor_r_region temp_AcceptorCorrect_r ];
            end
            
            subplot('position',[0.96 0.42 0.03 0.15]);
            x = -0.1:0.02:1.1;
            [hX,hN]=hist(TempFret12_region,x);
            barh(hN,hX,'k');
            temp=axis;
            temp(3)=-0.1;
            temp(4)=1.1;
            axis(temp);
            grid on;
            axis on;
            zoom on;
            
            subplot('position',[0.96 0.22 0.03 0.15]);
            x = -0.1:0.02:1.1;
            [hX,hN]=hist(TempFret13_region,x);
            barh(hN,hX,'k');
            temp=axis;
            temp(3)=-0.1;
            temp(4)=1.1;
            axis(temp);
            grid on;
            axis on;
            zoom on;
            
            subplot('position',[0.96 0.02 0.03 0.15]);
            x = -0.1:0.02:1.1;
            [hX,hN]=hist(TempFret23_region,x);
            barh(hN,hX,'k');
            temp=axis;
            temp(3)=-0.1;
            temp(4)=1.1;
            axis(temp);
            grid on;
            axis on;
            zoom on;
            
            i = i - 1;
        end
        
        if answer=='s'  %save region data
            again=1;
            filename_add_region = sprintf('%s_trace_%i.dat', filename_head, i); %XXX
            disp(filename_add_region);
            fid = fopen(filename_add_region, 'w');
            fprintf(fid, ' dbackground_g d2background_g abackground_g dbackground_r d2background_r abackground_r leakage12 leakage21 leakage13 leakage23 gamma12 gamma13 direct \n');
            fprintf(fid, '%f %f %f %f %f  %f %f %f %f %f  %f %f %f \n', dbackground_g, d2background_g, abackground_g, dbackground_r, d2background_r, abackground_r, leakage12, leakage21, leakage13, leakage23, gamma12, gamma13, direct);
            fprintf(fid, ' number firstpoint lastpoint length \n');
            output = double([ Temp_i; Tempfirstpoint; Templastpoint; Templength]);
            fprintf(fid, '%f %f %f %f \n', output);
            fclose(fid);
            
            output = [ Temptime_region' TempFret12_region' TempFret13_region' TempFret23_region' TempDonor_g_region' TempDonor2_g_region' TempAcceptor_g_region' TempDonor_r_region' TempDonor2_r_region' TempAcceptor_r_region'];
            save(filename_add_region_data,'output','-ascii');
        end
        
        disp('again end');
    end
    
    if answer=='l'
        [Xc,~,buttonc] = ginput(2); %ginput means graphicalinput from a mouse or cursor
        if (Xc(1)>Xc(2))
            temp = Xc(1);
            Xc(1) = Xc(2);
            Xc(2) = temp;
        end
        firstpoint = round(Xc(1)/(3*timeunit)); % fp means the Akaike final prediction for estimate
        lastpoint = round(Xc(2)/(3*timeunit)); % get database column privilige
        disp(firstpoint);
        disp(lastpoint);
        Temptime = [Temptime time_g(firstpoint:lastpoint)];
        TempFret12 = [TempFret12 Fret12(firstpoint:lastpoint)];
        TempFret13 = [TempFret13 Fret13(firstpoint:lastpoint)];
        TempFret23 = [TempFret23 Fret23(firstpoint:lastpoint)];
        TempDonor_g = [TempDonor_g DonorCorrect_g(firstpoint:lastpoint)];
        TempDonor2_g = [TempDonor2_g Donor2Correct_g(firstpoint:lastpoint)];
        TempDonor_r = [TempDonor_r DonorCorrect_r(firstpoint:lastpoint)];
        TempDonor2_r = [TempDonor2_r Donor2Correct_r(firstpoint:lastpoint)];
        TempAcceptor_g = [TempAcceptor_g AcceptorCorrect_g(firstpoint:lastpoint)];
        TempAcceptor_r = [TempAcceptor_r AcceptorCorrect_r(firstpoint:lastpoint)];
        output = [ Temptime' TempFret12' TempFret13' TempFret23' TempDonor_g' TempDonor2_g' TempAcceptor_g' TempDonor_r' TempDonor2_r' TempAcceptor_r'];
        file_outname = sprintf('%s_trace_%i.dat', filename_head, i);
        save(file_outname,'output','-ascii');
    end
    
    if answer=='l'
        [Xc,Yc,buttonc] = ginput(2); %ginput means graphicalinput from a mouse or cursor
        if (Xc(1)>Xc(2))
            temp = Xc(1);
            Xc(1) = Xc(2);
            Xc(2) = temp;
        end
        firstpoint = round(Xc(1)/(3*timeunit)); % fp means the Akaike final prediction for estimate
        lastpoint = round(Xc(2)/(3*timeunit)); % get database column privilige
        disp(firstpoint);
        disp(lastpoint);
        Temptime = [Temptime time_g(firstpoint:lastpoint)];
        TempFret12 = [TempFret12 Fret12(firstpoint:lastpoint)];
        TempFret13 = [TempFret13 Fret13(firstpoint:lastpoint)];
        TempFret23 = [TempFret23 Fret23(firstpoint:lastpoint)];
        TempDonor_g = [TempDonor_g DonorCorrect_g(firstpoint:lastpoint)];
        TempDonor2_g = [TempDonor2_g Donor2Correct_g(firstpoint:lastpoint)];
        TempDonor_r = [TempDonor_r DonorCorrect_r(firstpoint:lastpoint)];
        TempDonor2_r = [TempDonor2_r Donor2Correct_r(firstpoint:lastpoint)];
        TempAcceptor_g = [TempAcceptor_g AcceptorCorrect_g(firstpoint:lastpoint)];
        TempAcceptor_r = [TempAcceptor_r AcceptorCorrect_r(firstpoint:lastpoint)];
        output = [ Temptime' TempFret12' TempFret13' TempFret23' TempDonor_g' TempDonor2_g' TempAcceptor_g' TempDonor_r' TempDonor2_r' TempAcceptor_r'];
        file_outname = sprintf('%s_trace_%i.dat', filename_head, i);
        save(file_outname,'output','-ascii');
    end
    
    % calculate gamma factor (ashlee)
    if answer=='r'
        [raw_x, y] = ginput(4);
        x = round(raw_x / (3 * timeunit));
        % average intensities before change point
        
        d1d1_bef = mean(DonorCorrect_b(x(1):x(2)));
        d1d2_bef = mean(Donor2Correct_b(x(1):x(2)));
        d1ac_bef = mean(AcceptorCorrect_b(x(1):x(2)));
        
        d2d2_bef = mean(Donor2Correct_g(x(1):x(2)));
        d2ac_bef = mean(AcceptorCorrect_g(x(1):x(2)));
        
        d1d1_aft = mean(DonorCorrect_b(x(3):x(4)));
        d1d2_aft = mean(Donor2Correct_b(x(3):x(4)));
        d1ac_aft = mean(AcceptorCorrect_b(x(3):x(4)));
        
        d2d2_aft = mean(Donor2Correct_g(x(3):x(4)));
        d2ac_aft = mean(AcceptorCorrect_g(x(3):x(4)));
        
        delta_d1d1 = d1d1_aft - d1d1_bef;
        delta_d1d2 = d1d2_aft - d1d2_bef;
        delta_d1ac = d1ac_aft - d1ac_bef;
        
        delta_d2d2 = d2d2_aft - d2d2_bef;
        delta_d2ac = d2ac_aft - d2ac_bef;
        
        r12 = (-1) * delta_d1d1 / delta_d1d2;
        r13 = (-1) * delta_d1d1 / delta_d1ac;
        r23 = (-1) * delta_d2d2 / delta_d2ac;
        
        r12_list = [r12_list r12];
        r13_list = [r13_list r13];
        r23_list = [r23_list r23];
    end
    
    % calculate leakage (ashlee)
    
    if answer=='k'
        [raw_x, y] = ginput(2);
        x = round(raw_x / (3 * timeunit));
        d1d1 = mean(DonorCorrect_b(x(1):x(2)));
        d1d2 = mean(Donor2Correct_b(x(1):x(2)));
        d1ac = mean(AcceptorCorrect_b(x(1):x(2)));
        
        d2d2 = mean(Donor2Correct_g(x(1):x(2)));
        d2ac = mean(AcceptorCorrect_g(x(1):x(2)));
        
        l12 = d1d2 / d1d1;
        l13 = d1ac / d1d1;
        l23 = d2ac / d2d2;
        
        l12_list = [l12_list l12];
        l13_list = [l13_list l13];
        l23_list = [l23_list l23];
        
        fprintf('Leakage for this trace: l12 = %.2f, l13 = %.2f, l23 = %.2f', l12, l13, l23);
    end
    
    % calculate red laser direct excitation of cy7 (ashlee)
    if answer=='i'
        [raw_x, y] = ginput(4);
        x = round(raw_x / (3 * timeunit));
        % average intensities before change point
        d2d2_bef = mean(Donor2Correct_g(x(1):x(2)));
        d2ac_bef = mean(AcceptorCorrect_g(x(1):x(2)));
        
        d2d2_aft = mean(Donor2Correct_g(x(3):x(4)));
        d2ac_aft = mean(AcceptorCorrect_g(x(3):x(4)));
        
        dd2ac = (d2d2_aft + d2ac_aft) / (d2d2_bef + d2ac_bef);
        
        dd2ac_list = [dd2ac_list dd2ac];
    end
    
    % photobleaching analysis
    if answer == 'q'
        [x, y] = ginput(1);
        disp(x);
        t_list(i) = x;
    end
    
    % dwell time analysis
    if answer == 'w'
        disp('Click for beginning and end of states.');disp('Left/middle/right click for different states.');
        [time,y,button]=ginput;
        
        time1=time(button==1);
        for c=1:2:(sum(button==1)-1)
            t1=ceil(time1(c)/Timeunit);
            t2=ceil(time1(c+1)/Timeunit);
            DT1(end+1)=abs(time1(c+1)-time1(c));
            DT1a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
            DT1d(end+1)=mean(Donors(TracesCounter,t1:t2));
            DT1f(end+1)=mean(FRET_Time_Series(t1:t2));
        end
        time2=time(button==2);
        for c=1:2:sum(button==2)-1
            t1=ceil(time2(c)/Timeunit);t2=ceil(time2(c+1)/Timeunit);
            DT2(end+1)=abs(time2(c+1)-time2(c));
            DT2a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
            DT2d(end+1)=mean(Donors(TracesCounter,t1:t2));
            DT2f(end+1)=mean(FRET_Time_Series(t1:t2));
        end
        time3=time(button==3);
        for c=1:2:sum(button==3)-1
            t1=ceil(time3(c)/Timeunit);t2=ceil(time3(c+1)/Timeunit);
            DT3(end+1)=abs(time3(c+1)-time3(c));
            DT3a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
            DT3d(end+1)=mean(Donors(TracesCounter,t1:t2));
            DT3f(end+1)=mean(FRET_Time_Series(t1:t2));
        end
        
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
end

%%save region datas
if size(Temp_i)~=0
    fid = fopen(filename_add_region, 'w');
    fprintf(fid, ' dbackground_g d2background_g abackground_g dbackground_r d2background_r abackground_r leakage12 leakage21 leakage13 leakage23 gamma12 gamma13 direct \n');
    fprintf(fid, '%f %f %f %f %f  %f %f %f %f %f  %f %f %f \n', dbackground_g, d2background_g, abackground_g, dbackground_r, d2background_r, abackground_r, leakage12, leakage21, leakage13, leakage23, gamma12, gamma13, direct);
    fprintf(fid, ' number firstpoint lastpoint length \n');
    output = double([ Temp_i; Tempfirstpoint; Templastpoint; Templength]);
    fprintf(fid, '%f %f %f %f \n', output);
    fclose(fid);
    
    output = [ Temptime_region' TempFret12_region' TempFret13_region' TempFret23_region' TempDonor_g_region' TempDonor2_g_region' TempAcceptor_g_region' TempDonor_r_region' TempDonor2_r_region' TempAcceptor_r_region'];
    save(filename_add_region_data,'output','-ascii');
end
%%

disp('Program end.')

cd(OriginalDirectory);

if size(l12_list) ~= 0
    csvwrite('l12.csv', l12_list);
    csvwrite('l13.csv', l13_list);
    csvwrite('l23.csv', l23_list);
end

if size(r12_list) ~= 0
    csvwrite('r12.csv', r12_list);
    csvwrite('r13.csv', r13_list);
    csvwrite('r23.csv', r23_list);
end

if size(dd2ac_list) ~= 0
    csvwrite('direct_d2ac.csv', dd2ac_list);
end

%figure;
%photob_hist = histogram(t_list);
%xlabel('Time');
%ylabel('Count');
%saveas(photob_hist, [filename_head '_photobleaching_curve.png']);
%csvwrite([filename_head '_photobleaching.csv'], t_list);

% save dwell time data

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

clear all;
close all;
end
