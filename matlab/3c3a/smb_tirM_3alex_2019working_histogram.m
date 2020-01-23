% Single Molecule Biophysics Lab. in Seoul National University // MJ 2019 July
% edited by X. Feng Sep 4, 2019

%% path and filename setting
WorkingDirectory = pwd;

filenames = {'hel31' 'hel32' 'hel33' 'hel34' 'hel35'};
fret12 = []; % cy3-cy5 fret
intens1 = [];
intens2 = [];

for f = 1:5
    filename_head = char(filenames(f));
    folder_prefix = '';
    
    %% Correction parameters for FRET%%
    dbackground_b=0;
    d2background_b=0;abackground_b=0;
    
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
    
    direct = 0.1578; %ashlee: 0.1578;
    
    %% Options
    LaserOrderChange = 'y'; %Check this part when excitation laser order is matched.
    ColorNumber = 3;
    DyeType = 'cy235';      %not done: b->g, g->r, r->far red (color change only)
    DoseBinningNeed = 'n';  %#ok<*NASGU> % %binning required??
    binwidth=5;
    DoesFilterNeed = 'n';   % %filter required??
    DoseMovieNeed = 'n';    % Movie required??
    Is_Avg_and_E_save ='n'; % Average and E level save??
    FirstNumber = 10;       % histogram options
    LastNumber = 10;        % histogram options
    
    
    %% log file loading (time unit etc. Ha lab ver.)
    
    fileinfo = dir([filename_head '.log']);
    if sum(size(fileinfo)) == 1
        disp(['No log file : '  filename_head '.log']);
    end
    %date = fileinfo(1).date;
    fileid_log = fopen([filename_head '.log'],'r');		%% .log file
    %fileid_log = fopen(['hel1.log'],'r');
    logarray = textscan(fileid_log, '%s', 'Delimiter','\n');
    timeunit = 0.001*str2double(logarray{1,1}{strmatch('Exposure Time [ms]', logarray{1,1})+1});
    gain = str2double(logarray{1,1}{strmatch('Gain', logarray{1,1})+1});
    scaler = str2double(logarray{1,1}{strmatch('Data Scaler', logarray{1,1})+1});
    background_donor = str2double(logarray{1,1}{strmatch('Background', logarray{1,1})+1}); %#ok<*MATCH2>
    background_acceptor = str2double(logarray{1,1}{strmatch('Background', logarray{1,1})+1});
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
    fprintf('The length of the time traces is: %d\n', time_length);
    
    if mod(time_length,3)==0
        time_b = (0:+3:(time_length-3))*timeunit;
        time_g = (1:+3:(time_length-2))*timeunit;
        time_r = (2:+3:(time_length-1))*timeunit;
    end
    
    if mod(time_length,3)==1
        time_b = (0:+3:(time_length-3))*timeunit;
        time_g = (1:+3:(time_length-2))*timeunit;
        time_r = (2:+3:(time_length-1))*timeunit;
    end
    
    if mod(time_length,3)==2
        time_b = (0:+3:(time_length-3))*timeunit;
        time_g = (1:+3:(time_length-2))*timeunit;
        time_r = (2:+3:(time_length-1))*timeunit;
    end
    
    NumberofTraces = fread(fileid, 1, 'int16');
    NumberofPeaks = NumberofTraces/ColorNumber;
    fprintf('The number of traces and peaks are: %d, %d\n', NumberofTraces, NumberofPeaks);
    
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
    
    dbackground_b_temp = dbackground_b;
    d2background_b_temp = d2background_b;
    abackground_b_temp = abackground_b;
    
    dbackground_g_temp = dbackground_g;
    d2background_g_temp = d2background_g;
    abackground_g_temp = abackground_g;
    
    dbackground_r_temp = dbackground_r;
    d2background_r_temp = d2background_r;
    abackground_r_temp = abackground_r;
    
    
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
    % MJ edited
    
    DonorRawData_b = DonorRawData_1;
    Donor2RawData_b = Donor2RawData_1;
    AcceptorRawData_b = AcceptorRawData_1;
    DonorRawData_g = DonorRawData_2;
    Donor2RawData_g = Donor2RawData_2;
    AcceptorRawData_g = AcceptorRawData_2;
    DonorRawData_r = DonorRawData_0;
    Donor2RawData_r = Donor2RawData_0;
    AcceptorRawData_r = AcceptorRawData_0;
    
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
    
    if strcmp(DyeType, 'cy235') == 1
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
    
    
    DonorFirstData_g = reshape(DonorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    Donor2FirstData_g = reshape(Donor2RawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    AcceptorFirstData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    DonorLastData_g = reshape(DonorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
    Donor2LastData_g = reshape(Donor2RawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
    AcceptorLastData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
    DonorFirstData_r = reshape(DonorRawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    Donor2FirstData_r = reshape(Donor2RawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    AcceptorFirstData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    DonorLastData_r = reshape(DonorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
    Donor2LastData_r = reshape(Donor2RawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
    AcceptorLastData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
    
    
    DonorFirstData_g_after = reshape(DonorRawData_g(1:NumberofPeaks, 2:(FirstNumber+1)), NumberofPeaks*FirstNumber, 1);
    DonorFirstData_g_before = reshape(DonorRawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    Donor2FirstData_g_after = reshape(Donor2RawData_g(1:NumberofPeaks, 2:(FirstNumber+1)), NumberofPeaks*FirstNumber, 1);
    Donor2FirstData_g_before = reshape(Donor2RawData_g(1:NumberofPeaks, 1:FirstNumber), NumberofPeaks*FirstNumber, 1);
    
    E_level_12_after = Donor2FirstData_g_after./(Donor2FirstData_g_after + DonorFirstData_g_after);
    E_level_12_before = Donor2FirstData_g_before./(Donor2FirstData_g_before + DonorFirstData_g_before);
    
    E_level_12_onestep = [E_level_12(end)' E_level_12(1:(end-1))'];
    
    
    %% Start to servey
    t_list = -1 * ones(NumberofPeaks, 1);
    DT1=[];DT2=[];DT3=[];
    DT1a=[];DT2a=[];DT3a=[];
    DT1d=[];DT2d=[];DT3d=[];
    DT1f=[];DT2f=[];DT3f=[];
    DT1id = []; DT2id = []; DT3id = [];
    FirstSelectX=[];
    FirstSelectY=[];
    LastSelectX=[];
    LastSelectY=[];
    
    i=0;
    prev_i = -1;
    history_n = 0;
    history = zeros(1000, 1, 'int16');
    firstpoint = 1;
    lastpoint = time_length_each;
    junk=zeros(NumberofPeaks, 1);
    
    
    
    % Display traces
    
    while i < NumberofPeaks
        
        i = i + 1;
        
        %     if prev_i ~= i
        %         corrected = false;
        %     end
        
        if ColorNumber == 3
            if strcmp(DyeType, 'cy235') == 1
                DonorCorrect_b = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_b(i,:) - dbackground_b_temp) - leakage21 * (Donor2RawData_g(i,:) - d2background_b_temp));
                Donor2Correct_b = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_b(i,:) - d2background_b_temp) - leakage12 * (DonorRawData_b(i,:) - dbackground_b_temp));
                AcceptorCorrect_b = gamma13 * ((AcceptorRawData_b(i,:) - abackground_b_temp) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_b(i,:) - dbackground_b_temp) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_b(i,:) -d2background_b_temp));
                EachTotalCorrect_b = DonorCorrect_b + Donor2Correct_b + AcceptorCorrect_b;
                EachTotalCorrect_b = (EachTotalCorrect_b~=0).*EachTotalCorrect_b + (EachTotalCorrect_b==0)*1;	% remove zeros
                
                DonorCorrect_g = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_g(i,:) - dbackground_g_temp) - leakage21 * (Donor2RawData_g(i,:) - d2background_g_temp));
                Donor2Correct_g = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_g(i,:) - d2background_g_temp) - leakage12 * (DonorRawData_g(i,:) - dbackground_g_temp));
                AcceptorCorrect_g = gamma13 * ((AcceptorRawData_g(i,:) - abackground_g_temp) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_g(i,:) - dbackground_g_temp) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_g(i,:) -d2background_g_temp));
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
                %corrected = true;
            end
        end
        
        for j=2:time_length_each
            if (Fret13(j) < -0.2)
                Fret13(j) = 0;
            elseif (Fret13(j) > 1.2)
                Fret13(j) = 1;
            end
            if(Fret12(j) < -0.2)
                Fret12(j) = 0;
            elseif(Fret12(j) > 1.2)
                Fret12(j) = 1;
            end
            if(Fret23(j) < -0.2)
                Fret23(j) = 0;
            elseif(Fret23(j) > 1.2)
                Fret23(j) = 1;
            end
        end
        
        prev_i = i;
        
        fret12 = [fret12 mean(Fret12(1:10))];
        intens1_all = DonorCorrect_b + Donor2Correct_b + AcceptorCorrect_b;
        intens1 = [intens1 mean(intens1_all(1:10))];
        intens2_all = Donor2Correct_g + AcceptorCorrect_g;
        intens2 = [intens2 mean(intens2_all(1:10))];
        
        dbackground_b_temp = dbackground_b;
        d2background_b_temp = d2background_b;
        abackground_b_temp = abackground_b;
        
        dbackground_g_temp = dbackground_g;
        d2background_g_temp = d2background_g;
        abackground_g_temp = abackground_g;
        
        dbackground_r_temp = dbackground_r;
        d2background_r_temp = d2background_r;
        abackground_r_temp = abackground_r;
    end
end

molec_with_cy5 = false(N_mol, 1);
N_mol = length(intens2);
for i = 1:N_mol
    
    % has cy5
    if intens2(i) > 120 && intens2(i) < 450
        molec_with_cy5(i) = true;
    end
   
end

f1 = figure;
histogram(fret12(molec_with_cy5), 'BinWidth', 0.025);
xlim([-0.2 1.2]);
%ylim([0 600]);
title('Cy3 Cy5 FRET');
xlabel('FRET Efficiency');
ylabel('Count');
set(gca,'FontSize',20);

disp('Program end.')
