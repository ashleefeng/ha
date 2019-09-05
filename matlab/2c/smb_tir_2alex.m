
function smb_tir_2alex
% Single Molecule Biophysics Lab. in Seoul National University

clear all;
close all;
disp('           ')
%% Directory and filename
OriginalDirectory=cd;
WorkingDirectory='Z:\Ashlee\180731_B2_488SWR1_3color\2_B2_3000X';
cd(WorkingDirectory);
%% How many dye kind??  Two channel 2, three channel 3 %%
ColorNumber=3;
%% Filenames
filename_head = 'hel3';
filename_traces = [filename_head '.' num2str(ColorNumber) 'color_2alex_traces'];
filename_movie = [filename_head '.' num2str(ColorNumber) 'color_2alex_movies'];
filename_time = [filename_head '_time.dat'];

filename_add = [ filename_head '_select.dat'];
filename_add_region = [ filename_head '_region.dat' ];
filename_add_region_data=[ filename_head '_region_data.dat' ];
%% Data Correction %%
dbackground_g=0;
d2background_g=0;
abackground_g=0;

dbackground_r=0;
d2background_r=0;
abackground_r=0;

leakage12=0;   %0.11
leakage21=0;
leakage13=0;   %0.013
leakage23=0;   %0.12

gamma12=1;  %1
gamma13=1;  %2.36

direct = 0;   %0.19
%% binning required?? %%
DoseBinningNeed = 'n';
binwidth=5;
%% filter required?? %%
DoesFilterNeed = 'n';
%% Movie required?? %%
DoseMovieNeed = 'y';
%% Average and E level save?? %%
Is_Avg_and_E_save ='n';
%% Time unit is 'ms'? %%
Time_unit_ms = 'n';

BottomLimit_g=-100;
%UpperLimit_g=2800;
UpperLimit_g=1000;
BottomLimit_r=-100;
UpperLimit_r=1000;
%% Define time unit %%
fileinfo = dir([filename_head '.log']);
if sum(size(fileinfo)) == 1
    disp(['No log file : '  filename_head '.log']);
end
date = fileinfo.date;
fileid_log = fopen([filename_head '.log'],'r');		%% .log file
timeunit = textscan(fileid_log, '%*[Exposure time:] %f', 1);
timeunit = timeunit{1,1};
if Time_unit_ms =='y'
	timeunit= timeunit/1000;
end
 timeunit=0.1;     %% You want manual? Then use this line.
textscan(fileid_log, '%*[Acquisition mode:  Full 512x512 1x1 256x256 2x2 Binning]', 1);
gain = textscan(fileid_log, '%*[Gain:] %d', 1);
gain = gain{1,1};
scaler = textscan(fileid_log, '%*[Data scaler:] %f', 1);
scaler = scaler{1,1};
background = textscan(fileid_log, '%*[Background subtraction:] %d', 1);
background = background{1,1};
% background_acceptor = textscan(fileid_log, '%d', 1);
% background_acceptor = background_acceptor{1,1}
fclose(fileid_log);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading data %%
fileid = fopen(filename_traces,'r');	%% .trace file is made from idl(.run ana_all)
if fileid == -1
    disp(['No data file  '  filename_head]);
end
time_length = fread(fileid, 1, 'int32');
disp('The length of the time traces is: ')	% the length of the time trace
disp(time_length);
time = (0:(time_length-1))*timeunit;
if mod(time_length,2)==0
    time_g = (0:+2:(time_length-2))*timeunit;
    time_r = (1:+2:(time_length-1))*timeunit;
else
    time_g = (0:+2:(time_length-3))*timeunit;
    time_r = (1:+2:(time_length-2))*timeunit;
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
TempDonor_g_region = [];
TempDonor2_g_region = [];
TempAcceptor_g_region = [];
TempDonor_r_region = [];
TempDonor2_r_region = [];
TempAcceptor_r_region = [];

fileid_add_region = fopen(filename_add_region, 'r');
fileid_add_region_data = fopen(filename_add_region_data, 'r');
if fileid_add_region~=-1 && fileid_add_region_data~=-1
	fgetl(fileid_add_region);
	temp = fscanf(fileid_add_region, '%f %f %f %f %f  %f %f %f %f %f  %f %f %f\n', [13 1]);

	dbackground_g=double(temp(1,1));
    d2background_g=double(temp(2,1));
    abackground_g=double(temp(3,1));

    dbackground_r=double(temp(4,1));
    d2background_r=double(temp(5,1));
    abackground_r=double(temp(6,1));

    leakage12=double(temp(7,1));   %0.11
    leakage21=double(temp(8,1));
    leakage13=double(temp(9,1));   %0.013
    leakage23=double(temp(10,1));   %0.12

    gamma12=double(temp(11,1));  %1
    gamma13=double(temp(12,1));  %2.36

    direct=double(temp(13,1));   %0.19

	fgetl(fileid_add_region);
	temp = fscanf(fileid_add_region, '%f %f %f %f\n', [4 inf]);
	Temp_i = int32(temp(1,:));
	Tempfirstpoint = int32(temp(2,:));
	Templastpoint = int32(temp(3,:));
	Templength = int32(temp(4,:));

	temp = fscanf(fileid_add_region_data, '%f %f %f %f %f  %f %f %f %f %f\n', [10 inf]);
    Temptime_region = temp(1,:);
    TempFret12_region = temp(2,:);
    TempFret13_region = temp(3,:);
    TempFret23_region = temp(4,:);
    TempDonor_g_region = temp(5,:);
    TempDonor2_g_region = temp(6,:);
    TempAcceptor_g_region = temp(7,:);
    TempDonor_r_region = temp(8,:);
    TempDonor2_r_region = temp(9,:);
    TempAcceptor_r_region = temp(10,:);
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
time_length_each = floor(time_length/2);
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

% Rawdata를 나누어서 trace만들기
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
            DonorRawData_1(i,j) = Data(i*3-2,j*2-1);
            Donor2RawData_1(i,j) = Data(i*3-1,j*2-1);
            AcceptorRawData_1(i,j) = Data(i*3,j*2-1);
            DonorRawData_2(i,j) = Data(i*3-2,j*2);
            Donor2RawData_2(i,j) = Data(i*3-1,j*2);
            AcceptorRawData_2(i,j) = Data(i*3,j*2);
        end
    end
end

clear Data;



%% calculate, plot and save average traces %%
%각 trace의 average를 구한후, 첫번째 laser excitation의 Cy3 channel과 두번째 laser
%excitation의 cy3 channel을 비교한 후, 큰 쪽의 trace를 green laser excitation으로 작은 쪽의
%trace를 red laser excitation으로 배정
if sum(sum(DonorRawData_1, 1)) > sum(sum(DonorRawData_2, 1))
    DonorRawData_g = DonorRawData_1;
    Donor2RawData_g = Donor2RawData_1;
    AcceptorRawData_g = AcceptorRawData_1;
    DonorRawData_r = DonorRawData_2;
    Donor2RawData_r = Donor2RawData_2;
    AcceptorRawData_r = AcceptorRawData_2;
    color_order = 0;
else
    DonorRawData_g = DonorRawData_2;
    Donor2RawData_g = Donor2RawData_2;
    AcceptorRawData_g = AcceptorRawData_2;
    DonorRawData_r = DonorRawData_1;
    Donor2RawData_r = Donor2RawData_1;
    AcceptorRawData_r = AcceptorRawData_1;
	color_order = 1;
end

DonorRawAverage_g = sum(DonorRawData_g, 1) / NumberofPeaks;
Donor2RawAverage_g = sum(Donor2RawData_g, 1) /NumberofPeaks;
AcceptorRawAverage_g = sum(AcceptorRawData_g, 1) / NumberofPeaks;
DonorRawAverage_r = sum(DonorRawData_r, 1) / NumberofPeaks;
Donor2RawAverage_r = sum(Donor2RawData_r, 1) /NumberofPeaks;
AcceptorRawAverage_r = sum(AcceptorRawData_r, 1) / NumberofPeaks;

clear DonorRawData_1;
clear Donor2RawData_1;
clear AcceptorRawData_1;
clear DonorRawData_2;
clear Donor2RawData_2;
clear AcceptorRawData_2;

figure('Name','Raw Data Ensemble Average');
hdl1 = gcf;

%Green excitation의 signal average
subplot(2,1,1);
plot(time_g, DonorRawAverage_g - dbackground_g, 'g', time_g, Donor2RawAverage_g - d2background_g, 'r', time_g, AcceptorRawAverage_g - abackground_g, 'b');
title('Raw data Average donor, donor2 and acceptor signal in Green excitation');
zoom on;

%Red excitation의 signal average
subplot(2,1,2);
plot(time_r, DonorRawAverage_r - dbackground_r, 'g', time_r, Donor2RawAverage_r - d2background_r, 'r', time_r, AcceptorRawAverage_r - abackground_r, 'b');
title('Raw data Average donor, donor2 and acceptor signal in Red excitation');
zoom on;

if Is_Avg_and_E_save =='y'
    AverageOutput_g = [time_g' DonorRawAverage_g' Donor2RawAverage_g' AcceptorRawAverage_g'];
    AverageOutput_r = [time_r' DonorRawAverage_r' Donor2RawAverage_r' AcceptorRawAverage_r'];
    save([filename_head '_avg_g.dat'], 'AverageOutput_g', '-ascii');
    save([filename_head '_avg_r.dat'], 'AverageOutput_r', '-ascii');
end

%% calculate E level from the first 10 points and plot histograms of E level and total intensity. Also save the same info
FirstStart = 11;
FirstEnd = 20;
LastNumber = 10;
tempDonor_g = reshape(DonorRawData_g(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
tempDonor2_g = reshape(Donor2RawData_g(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
tempAcceptor_g = reshape(AcceptorRawData_g(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
tempDonor_r = reshape(DonorRawData_r(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
tempDonor2_r = reshape(Donor2RawData_r(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
tempAcceptor_r = reshape(AcceptorRawData_r(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);

%Fret23 = AcceptorCorrect_r./(Donor2Correct_r+AcceptorCorrect_r);
%EachTotalCorrect_g = DonorCorrect_g + Donor2Correct_g + AcceptorCorrect_g;
%EachTotalCorrect_g = (EachTotalCorrect_g~=0).*EachTotalCorrect_g + (EachTotalCorrect_g==0)*1;	% remove zeros
%Fret12 = Donor2Correct_g./((1-Fret23).*DonorCorrect_g + Donor2Correct_g);
%Fret13 = (AcceptorCorrect_g - Fret23.*(Donor2Correct_g + AcceptorCorrect_g))./(DonorCorrect_g + AcceptorCorrect_g - Fret23 .* (EachTotalCorrect_g));

EachTotal_23_r = tempDonor2_r + tempAcceptor_r;
EachTotal_23_r = (EachTotal_23_r~=0).*EachTotal_23_r + (EachTotal_23_r==0)*1;	% remove zeros
EachTotal_123_g = tempDonor_g + tempDonor2_g + tempAcceptor_g;
EachTotal_123_g = (EachTotal_123_g~=0).*EachTotal_123_g + (EachTotal_123_g==0)*1;	% remove zeros

E_level_23 = tempAcceptor_r./EachTotal_23_r;
E_level_12 = tempDonor2_g./((1-E_level_23).*tempDonor_g + tempDonor_g);
E_level_13 = (tempAcceptor_g - E_level_23.*(tempDonor2_g + tempAcceptor_g))./(tempDonor_g + tempAcceptor_g - E_level_23.*(EachTotal_123_g));

if Is_Avg_and_E_save =='y'
    E_level_output_g = [E_level_g EachTotal_g];
    E_level_output_r = [E_level_r EachTotal_r];
    save([filename_head '_elevel_10p_g.dat'],'E_level_output_g','-ascii');
    save([filename_head '_elevel_10p_r.dat'],'E_level_output_r','-ascii');
end

figure('Name','Raw data analysis');
hdl2 = gcf;

% subplot(3,4,1); % 2*3 figure_ upper left (the last number shows the location of figure)
% hist(E_level_12,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
% temp=axis;
% temp(1)=-0.1;
% temp(2)=1.1;
% axis(temp);
% title([ 'first ' num2str(FirstEnd) 'p Raw FRET_12 histogram' ]);
% zoom on;

% subplot(3,4,2); % 2*3 figure_ upper left (the last number shows the location of figure)
% hist(E_level_13,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
% temp=axis;
% temp(1)=-0.1;
% temp(2)=1.1;
% axis(temp);
% title([ 'first ' num2str(FirstEnd) 'p Raw FRET_13 histogram' ]);
% zoom on;

subplot(3,4,3); % 2*3 figure_ upper left (the last number shows the location of figure)
hist(E_level_23,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
axis(temp);
title([ 'first ' num2str(FirstEnd) 'p Raw FRET_23 histogram' ]);
zoom on;

subplot(3,4,4);
hist(EachTotal_123_g,-100:50:4000);
temp=axis;
temp(1)=-100;
temp(2)=4000;
axis(temp);
title([ 'first ' num2str(FirstEnd) 'p Raw Total intensity histogram' ]);
zoom on;

% subplot(2,2,3);
% plot(E_level_12, EachTotal_123_g,'b+', 'MarkerSize', 2);
% temp=axis;
% temp(1)=-0.1;
% temp(2)=1.1;
% temp(3)=-100;
% temp(4)=4000;
% axis(temp);
% title([ 'first ' num2str(FirstEnd) 'p Raw Total Intensity_123 vs. FRET' ]);
% zoom on;

DonorFirstData_g = reshape(DonorRawData_g(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
Donor2FirstData_g = reshape(Donor2RawData_g(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
AcceptorFirstData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
size(DonorRawData_g)
time_length_each
DonorLastData_g = reshape(DonorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
Donor2LastData_g = reshape(Donor2RawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
AcceptorLastData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);

DonorFirstData_r = reshape(DonorRawData_r(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
Donor2FirstData_r = reshape(Donor2RawData_r(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
AcceptorFirstData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, FirstStart:FirstEnd), NumberofPeaks*(FirstEnd - FirstStart + 1), 1);
DonorLastData_r = reshape(DonorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
Donor2LastData_r = reshape(Donor2RawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
AcceptorLastData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);


subplot(5,4,11);
hist(DonorFirstData_g,-300:5:UpperLimit_g);
temp=axis;
temp(1)=-300;
temp(2)=UpperLimit_g;
axis(temp);
title([ 'first ' num2str((FirstEnd - FirstStart + 1)) 'p Raw Donor Intensity histogram']);
zoom on;

subplot(5,4,12);
hist(DonorLastData_g,-300:5:UpperLimit_g);
temp=axis;
temp(1)=-300;
temp(2)=UpperLimit_g;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Donor Intensity histogram']);
zoom on;

subplot(5,4,15);
hist(Donor2FirstData_g,-300:5:UpperLimit_g);
temp=axis;
temp(1)=-300;
temp(2)=UpperLimit_g;
axis(temp);
title([ 'first ' num2str((FirstEnd - FirstStart + 1)) 'p Raw Donor2 Intensity histogram']);
zoom on;

subplot(5,4,16);
hist(Donor2LastData_g,-300:5:UpperLimit_g);
temp=axis;
temp(1)=-300;
temp(2)=UpperLimit_g;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Donor2 Intensity histogram']);
zoom on;

subplot(5,4,19);
hist(AcceptorFirstData_g,-300:5:UpperLimit_g);
temp=axis;
temp(1)=-300;
temp(2)=UpperLimit_g;
axis(temp);
title([ 'first ' num2str((FirstEnd - FirstStart + 1)) 'p Raw Acceptor Intensity histogram']);
zoom on;

subplot(5,4,20);
hist(AcceptorLastData_g,-300:5:UpperLimit_g);
temp=axis;
temp(1)=-300;
temp(2)=UpperLimit_g;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Acceptor Intensity histogram']);
zoom on;



%% Start to servey
Temptime = [];
%TempFret = [];
%TempDonor = [];
%TempAcceptor = [];
TempFret12 = [];
TempFret13 = [];
TempFret23 = [];
TempDonor_g = [];
TempDonor2_g = [];
TempAcceptor_g = [];
TempDonor_r = [];
TempDonor2_r= [];
TempAcceptor_r = [];

figure('Name','Trace analysis');
hdl_trace=gcf;

i=0;
history_n = 0;
history = zeros(1000, 1, 'int16');
while i < NumberofPeaks
    i = int32(int32(i) + 1);

    again=1;
	while ((again==1) && (i <= NumberofPeaks))
        if ColorNumber == 3
            DonorCorrect_g = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_g(i,:) - dbackground_g) - leakage21 * (Donor2RawData_g(i,:) - d2background_g));
            Donor2Correct_g = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_g(i,:) - d2background_g) - leakage12 * (DonorRawData_g(i,:) - dbackground_g)); 
    		AcceptorCorrect_g = gamma13 * ((AcceptorRawData_g(i,:) - abackground_g) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_g(i,:) - dbackground_g) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_g(i,:) -d2background_g));
			EachTotalCorrect_g = DonorCorrect_g + Donor2Correct_g + AcceptorCorrect_g;
			EachTotalCorrect_g = (EachTotalCorrect_g~=0).*EachTotalCorrect_g + (EachTotalCorrect_g==0)*1;	% remove zeros
			
			DonorCorrect_r = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_r(i,:) - dbackground_r) - leakage21 * (Donor2RawData_r(i,:) - d2background_r));
			Donor2Correct_r = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_r(i,:) - d2background_r) - leakage12 * (DonorRawData_r(i,:) - dbackground_r)); 
			AcceptorCorrect_r = gamma13 * ((AcceptorRawData_r(i,:) - abackground_r) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_r(i,:) - dbackground_r) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_r(i,:) -d2background_r));
			AcceptorCorrect_r = AcceptorCorrect_r - direct * (Donor2Correct_r + AcceptorCorrect_r);
			EachTotalCorrect_r = Donor2Correct_r + AcceptorCorrect_r;
			EachTotalCorrect_r = (EachTotalCorrect_r~=0).*EachTotalCorrect_r + (EachTotalCorrect_r==0)*1;	% remove zeros

			Fret23 = AcceptorCorrect_r./EachTotalCorrect_r;
			Fret12 = Donor2Correct_g./((1-Fret23).*DonorCorrect_g + Donor2Correct_g);
			Fret13 = (AcceptorCorrect_g - Fret23.*(Donor2Correct_g + AcceptorCorrect_g))./(DonorCorrect_g + AcceptorCorrect_g - Fret23 .* (EachTotalCorrect_g));
        end
        if ColorNumber == 2
            DonorCorrect_g = (DonorRawData_g(i,:) - dbackground_g) + leakage12 * (DonorRawData_g(i,:) - dbackground_g);
            Donor2Correct_g = gamma12 * ((Donor2RawData_g(i,:) - d2background_g) - leakage12 * (DonorRawData_g(i,:) - dbackground_g));
            AcceptorCorrect_g = 0 * DonorRawData_g(i,:);
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
	figure(hdl_trace);
	%subplot('position',[0.1 0.82 0.8 0.15]); %first trace : Green excitation 3-color corrected trace
	subplot(3, 1, 1);
    plot(time_g, DonorCorrect_g, 'b', time_g, Donor2Correct_g, 'g', time_g, AcceptorCorrect_g, 'r');
	temp=axis;
	temp(3)=BottomLimit_g;
	%temp(4)=UpperLimit_g;
	grid on;
	axis(temp);
	%title(['                    Molecule ' num2str(i) '  / ' num2str(NumberofPeaks) '              file: ' filename_head ]); %'   cor.: ' num2str(dbackground) '  ' num2str(abackground) '  ' num2str(d2background) ' ' num2str(leakage12) '  ' num2str(leakage13) ' ' num2str(leakage23) ' ' num2str(gamma12) ' ' num2str(gamma13) '   filter: ' DoesFilterNeed ' Fret:' num2str(FretEc_filter) ' Mean:' num2str(MeanRange) ' Acceptor:' num2str(Acceptor_filter) ' Total:' num2str(Total_filter) ]);
	title(['Blue Laser Molecule  ' num2str(i) '  / ' num2str(NumberofPeaks) '       file:' filename_head  ' cor.: ' num2str(dbackground_g) '  ' num2str(d2background_g) '  ' num2str(abackground_g) '  ' num2str(dbackground_r) '  ' num2str(d2background_r) '  ' num2str(abackground_r) '  ' num2str(leakage12) '  ' num2str(leakage21) '  ' num2str(leakage13) '  ' num2str(leakage23) '  ' num2str(gamma12) '  ' num2str(gamma13) ]);
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
	%subplot('position',[0.1 0.62 0.8 0.15]);
    subplot(3, 1, 2);
    plot(time_g, DonorCorrect_r, 'b', time_g, Donor2Correct_r, 'g', time_g, AcceptorCorrect_r, 'r');  
	temp=axis;
	temp(3)=BottomLimit_r;
	%temp(4)=UpperLimit_r;
	grid on;
	axis(temp);
	title(['Green Laser Molecule  ' num2str(i) '  / ' num2str(NumberofPeaks) '            timeunit: ' num2str(timeunit) '  gain: ' num2str(gain) '  scaler: ' num2str(scaler)  '  spot diameter: ' num2str(SpotDiameter) '    ' date]);
    zoom on;

	
% 	subplot('position',[0.93 0.42 0.03 0.15]);
% 	x = -0.1:0.02:1.1;
% 	[hX,hN]=hist(Fret12,x);
%     barh(hN,hX,'k');
% 	temp=axis;
% 	temp(3)=-0.1;
% 	temp(4)=1.1;
% 	axis(temp);
% 	grid on;
% 	axis on;
% 	zoom on;
	
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
    
% 	subplot('position',[0.1 0.42 0.8 0.15]);
% %	FretEc=(1./(1+gamma*(donorcorrect(i,:)./acceptorcorrect(i,:))));
% 	hFretLine = plot(time_g, Fret12, FirstSelectX, FirstSelectY, LastSelectX, LastSelectY, bintime, binEcorrect, 'k');
% 	temp=axis;
% 	temp(3)=-0.1;
% 	temp(4)=1.1;
% 	axis(temp);
% 	grid on;
% 	zoom on;

% 	subplot('position',[0.1 0.22 0.8 0.15]);
%     subplot(3, 1, 3);
% 	plot(time_g, Fret13, bintime, binEraw, 'k');
% 	temp=axis;
% 	temp(3)=-0.1;
% 	temp(4)=1.1;
% 	axis(temp);
% 	grid on;
% 	zoom on;
    
% 	subplot('position',[0.93 0.22 0.03 0.15]);
% 	x = -0.1:0.02:1.1;
% 	[hX,hN]=hist(Fret13,x);
% 	barh(hN,hX,'k');
% 	temp=axis;
% 	temp(3)=-0.1;
% 	temp(4)=1.1;
% 	axis(temp);
% 	grid on;
% 	axis on;
% 	zoom on;

    subplot(3, 1, 3);
    plot(time_g, Fret23, bintime, binEraw, 'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp)
    grid on;
    zoom on;
% 
% 	subplot('position',[0.93 0.02 0.03 0.15]);
% 	x = -0.1:0.02:1.1;
% 	[hX,hN]=hist(Fret23,x);
% 	barh(hN,hX,'k');
% 	temp=axis;
% 	temp(3)=-0.1;
% 	temp(4)=1.1;
% 	axis(temp);
% 	grid on;
% 	axis on;
% 	zoom on;    
    
	if DoseMovieNeed == 'y'
		fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
		%startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + color_order*peak_height*peaks_total_width;
        Xpoint = 6;
        startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + color_order*peak_height*peaks_total_width;
		for j=1:peak_height
			fseek(fileid_movie, startpoint , 'bof');
			peak_line=fread(fileid_movie, peak_height*3, 'uint8');
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
        
		fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
		%startpoint = uint32(4 + (i-1)*ColorNumber*peak_height) + (1-color_order)*peak_height*peaks_total_width;
        startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + (1-color_order)*peak_height*peaks_total_width;
		for j=1:peak_height
			fseek(fileid_movie, startpoint , 'bof');
            peak_line=fread(fileid_movie, peak_height*3, 'uint8');
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
	end


	again=1;
	while again==1
		again=0;
		disp([num2str(i) ' (n= next region, l=save select, u=select region, d=delete region, s=save region, h=histogram select, p=display movie)']);
		keyanswer =input('(t=terminate program, b=back, g=go) : ','s');
		answer = sscanf(keyanswer, '%s %*s');
		numberofanswer = sscanf(keyanswer, '%*s %f');

		
		if answer=='n'
			sum(find( i < Temp_i, 1, 'first' ))
			if sum(find( i < Temp_i, 1, 'first' )) ~= 0
				i = Temp_i(find( i < Temp_i, 1, 'first' ))-1;
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
				fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
				
				startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + color_order*peak_height*peaks_total_width;
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
                title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
        
				fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + (1-color_order)*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*3, 'uint8');
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
		[Xc,Yc,buttonc] = ginput(2); %ginput means graphicalinput from a mouse or cursor
		if (Xc(1)>Xc(2))
			temp = Xc(1);
			Xc(1) = Xc(2);
			Xc(2) = temp;
		end
		firstpoint = round(Xc(1)/(2*timeunit)); % fp means the Akaike final prediction for estimate
		lastpoint = round(Xc(2)/(2*timeunit)); % get database column privilige
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
		save(filename_add,'output','-ascii');
	end

    if answer=='s'
        filename_add = [ filename_head '_trace.dat'];
        output = [ time_g' DonorCorrect_g' Donor2Correct_g' AcceptorCorrect_g' DonorCorrect_r' Donor2Correct_r' AcceptorCorrect_r' Fret12' Fret13' Fret23'];
        save(filename_add, 'output', '-ascii');
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

clear all;
close all;
end
