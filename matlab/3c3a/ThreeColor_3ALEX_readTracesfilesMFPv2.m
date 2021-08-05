
close all;
fclose('all');

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

leakage12=0.1066;   %0.11
leakage21=0.0;
leakage13=0.0083;   %0.013
leakage23=0.0446;   %0.12

gamma12=0.8730;  %1
gamma23 = 2.62; 
gamma13=gamma12*gamma23;

direct = 0;


%%%read data 
pth=input('directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
    pth=pwd;
end
cd(pth);
save_file=pth;

disp(pth);
FolderDir=dir;
[nf,dumn]=size(FolderDir);

fname=input('index # of filename [default=1]  ');
if isempty(fname)
   fname=1;
end
fname=num2str(fname);
disp(['hel' fname '.traces']);
disp(['film' fname '.traces']);

timeunit=input('time unit [default=0.05 sec]  ');
if isempty(timeunit)'blah'
    timeunit=0.2;
end

%select the folder to which the files need to be saved
newfolder = ['selected traces'];
mkdir(newfolder);

fid=fopen(['hel' fname '_raw' '.traces'],'r');
%fid=fopen(['film' fname '.traces'],'r');

fid2=fopen(['hel' fname '_bsl' '.traces'],'r');

%first line of binary file specifies length of trace
len=fread(fid,1,'int32');
disp('The len of the time traces is: ')
disp(len);

%first line of binary file specifies length of trace
len=fread(fid2,1,'int32');
disp('The len of the time traces2 is: ')
disp(len);

%number of traces
Ntraces=fread(fid,1,'int16');
disp('The number of traces is: ')
disp(Ntraces/2);

%number of traces
Ntraces=fread(fid2,1,'int16');
disp('The number of traces2 is: ')
disp(Ntraces/2);

%raw is a linear array, looks like it
raw=fread(fid,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid);

%raw is a linear array, looks like it
raw2=fread(fid2,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid2);

%convert background into traces
index2=(1:Ntraces*len);
Data2=zeros(Ntraces,len);
intensity32=zeros(Ntraces/3,len);
intensity52=zeros(Ntraces/3,len);
intensity72=zeros(Ntraces/3,len);
% total=zeros(Ntraces/3,1);
Data2(index2)=raw2(index2);
for i=1:(Ntraces/3)
    intensity32(i,:)=Data2(i*3-2,:);
    intensity52(i,:)=Data2(i*3-1,:);
    intensity72(i,:)=Data2(i*3,:);
%     tempD=sum(Cy3(i,(11:20)),2);
%     tempA=sum(Cy5(i,(11:20)),2);
%     total(i)=(tempA+tempD)/10.; 
end

%convert data into traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
intensity3=zeros(Ntraces/3,len);
intensity5=zeros(Ntraces/3,len);
intensity7=zeros(Ntraces/3,len);
% total=zeros(Ntraces/3,1);
Data(index)=raw(index);
for i=1:(Ntraces/3)
    intensity3(i,:)=Data(i*3-2,:);
    intensity5(i,:)=Data(i*3-1,:);
    intensity7(i,:)=Data(i*3,:);
%     tempD=sum(Cy3(i,(11:20)),2);
%     tempA=sum(Cy5(i,(11:20)),2);
%     total(i)=(tempA+tempD)/10.; 
end

Data_b= Data-Data2;
Data_b(:,501)=0;  %adjusted for datasets not divisible by 3
 
num_frames = size(Data_b,2);


for i = 0:(num_frames-1)/3
    k=i+1;
greenlaser_time(k)= (i*3+1)*timeunit;
redlaser_time(k) = (i*3+2)*timeunit;
farredlaser_time(k) = (i*3+3)*timeunit;
    end

    
    

% figure;
% hist(total,80);
% grid on;
% zoom on;
% 
% fcutoff1=input('low cutoff intensity: ','s');
% cutoff1=str2double(fcutoff1);
% if isempty(cutoff1)
%     cutoff1=300;
% end
% fcutoff2=input('high cutoff intensity: ','s');
% cutoff2=str2double(fcutoff2);
% if isempty(cutoff2)
%     cutoff2=2500;
% end

% i=0;
% n=0;
 N_mol=Ntraces/3;

 len=len+1; %adjusted for number of frames not divisible by 3
 
len_each = len/3;

 for i=1:N_mol
        for j=1:len_each
            DonorRawData_0(i,j) = Data_b(i*3-2,j*3-2);
            Donor2RawData_0(i,j) = Data_b(i*3-1,j*3-2);
            AcceptorRawData_0(i,j) = Data_b(i*3,j*3-2);
            DonorRawData_1(i,j) = Data_b(i*3-2,j*3-1);
            Donor2RawData_1(i,j) = Data_b(i*3-1,j*3-1);
            AcceptorRawData_1(i,j) = Data_b(i*3,j*3-1);
            DonorRawData_2(i,j) = Data_b(i*3-2,j*3);
            Donor2RawData_2(i,j) = Data_b(i*3-1,j*3);
            AcceptorRawData_2(i,j) = Data_b(i*3,j*3);
        end
 end
    
 %% 
 % Corrections

Donor2Correct_0 = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_0 - d2background_0) - leakage12 * (DonorRawData_0 - dbackground_0));
DonorCorrect_0 = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * (DonorRawData_0 - dbackground_0) - leakage21 * (Donor2RawData_0 - d2background_0);
AcceptorCorrect_0 = gamma13 * ((AcceptorRawData_0 - abackground_0) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_0 - dbackground_0) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_0 - d2background_0));
EachTotalCorrect_0 = DonorCorrect_0 + Donor2Correct_0 + AcceptorCorrect_0;
EachTotalCorrect_0 = (EachTotalCorrect_0~=0).*EachTotalCorrect_0 + (EachTotalCorrect_0==0)*1;	% remove zeros

DonorCorrect_1 = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * (DonorRawData_1 - dbackground_1) - leakage21 * (DonorRawData_1 - d2background_1);
Donor2Correct_1 = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * (Donor2RawData_1 - d2background_1) - leakage12 * (DonorRawData_1 - dbackground_1);
AcceptorCorrect_1 = gamma13 * (AcceptorRawData_1 - abackground_1) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_1 - dbackground_1) + ((leakage13 * leakage21 - leakage23) / (1 - leakage12 * leakage21)) * (Donor2RawData_1 -d2background_1);
AcceptorCorrect_1 = AcceptorCorrect_1 - direct * (Donor2Correct_1 + AcceptorCorrect_1);
EachTotalCorrect_1 = Donor2Correct_1 + AcceptorCorrect_1;

AcceptorCorrect_2 =  AcceptorRawData_2 - abackground_2;
DonorCorrect_2=DonorRawData_2;
Donor2Correct_2=Donor2RawData_2;

fret12 = AcceptorCorrect_1 ./ EachTotalCorrect_1;
fret01 = Donor2Correct_0 ./ ((1-fret12).*DonorCorrect_0 + Donor2Correct_0);
fret02 = (AcceptorCorrect_0 - fret12.*(Donor2Correct_0 + AcceptorCorrect_0)) ./ (DonorCorrect_0 + AcceptorCorrect_0 - fret12 .* (EachTotalCorrect_0));

global Donor2Correct_0 DonorCorrect_0 AcceptorCorrect_0 EachTotalCorrect_0 DonorCorrect_1... 
    Donor2Correct_1 AcceptorCorrect_1 EachTotalCorrect_1 AcceptorCorrect_2 DonorCorrect_2...
    Donor2Correct_2 fret12 fret01 fret02 greenlaser_time redlaser_time farredlaser_time i

%% 
% while (N_mol-i) > 0
%     i = i+1;
%     if total(i) < cutoff1 || total(i) > cutoff2
%         n = n+1;
%     end
% end
% 
%  L12=0.08;
%  L13=0.00;
%  L23=0.0325;

% %Frame(end-9) to Frame(end) are green laser for BG_GE3
% %Frame(end-19) to Frame(end-10) are red laser for BG_RE
% GreenChannel_BG_GE = mean(intensity3(:,end-9:end),2);
% RedChannel_BG_GE = mean(intensity5(:,end-9:end),2);
% FarRedChannel_BG_GE = mean(intensity7(:,end-9:end),2);
% GreenChannel_BG_RE = mean(intensity3(:,end-19:end-10),2);
% RedChannel_BG_RE = mean(intensity5(:,end-19:end-10),2);
% FarRedChannel_BG_RE = mean(intensity7(:,end-19:end-10),2);
% 
% intensity3_GE=intensity3(:,11:2:end-20)-median(GreenChannel_BG_GE);
% intensity5_GE=intensity5(:,11:2:end-20)-median(RedChannel_BG_GE);
% intensity7_GE=intensity7(:,11:2:end-20)-median(FarRedChannel_BG_GE);
% Time_GE=time(:,11:2:end-20);
% 
% intensity3_RE=intensity3(:,12:2:end-20)-min(GreenChannel_BG_RE);
% intensity5_RE=intensity5(:,12:2:end-20)-min(RedChannel_BG_RE);
% intensity7_RE=intensity7(:,12:2:end-20)-min(FarRedChannel_BG_RE);
% Time_RE=time(:,12:2:end-20);
% 
% intensity3_RE_2=intensity3(:,12:2:end-20)-min(GreenChannel_BG_RE);
% intensity5_RE_2=intensity5(:,12:2:end-20)-min(RedChannel_BG_RE);
% intensity7_RE_2=intensity7(:,12:2:end-20)-min(FarRedChannel_BG_RE);
% 
% choice=[];
% mkdir('Leakage_Background_Corrected');

hdl=figure;
i=0;
t=0;

counter_1=0;
counter_transition=0;


for i=1:N_mol
    
figure(hdl);
    GE=subplot(6,1,1);
    %A = get(gca,'position');          % gca points at the second one
     %set(gca,'unit','pixels')
    %A(1,4) = A(1,4) * 1.2;              % change the height 
    %A(1,2) = A(1,2) + .2*A(1,4);         % change the vertical position
    %set(gca,'position',A);            % set the values you just changed
    %set(gca,'unit','pixels')
    %GE.Position = [.13 .7 .83 .01];
    %GE.Position = [.13 .7 .83 .079];
    %GE.Position = [.13 .7 .775 .15] % define your position
    I1=(DonorCorrect_0(i,:));
    I2=((Donor2Correct_0(i,:)));
    I3=(AcceptorCorrect_0(i,:));
    plot(greenlaser_time,I1,'g',greenlaser_time,I2,'r', greenlaser_time,I3,'k', 'LineWidth', 1);
     ylim([-50 600]);
    %title(['  Molecule ' num2str(n) ' of ' num2str(N_mol) ' GE']);
    %xlabel('Time(s)');
    %ylabel('Intensity (A.U.)');

% 
% while n < N_mol
%     
%     if choice==0
%         n=n-1;
%     else
%         n=n+1;
%     end
    
%     if total(i) < cutoff1 || total(i) > cutoff2
%         continue;
%     end
    
    %trace window
   

    figure(hdl);
        %Afig=get(hdl,'position');
%Afig(1,4) = Afig(1,4) * 3;
%set(hdl,'position',Afig);
    RE=subplot(6,1,2);
    %A = get(gca,'position');          % gca points at the second one
    %set(gca,'unit','pixels')
    %A(1,4) = A(1,4) * 1.2;              % change the height 
    %A(1,2) = A(1,2) - .2*A(1,4);         % change the vertical position
    %set(gca,'position',A);            % set the values you just changed
    %RE.Position = [.13 .7 .775 .079]; % define your position
    I1=(DonorCorrect_1(i,:));
    I2=((Donor2Correct_1(i,:)));
    I3=(AcceptorCorrect_1(i,:));
    plot(redlaser_time,I1,'g',redlaser_time,I2,'r', redlaser_time,I3,'k');
     ylim([-50 800]);
    %title(['  Molecule ' num2str(n) ' of ' num2str(N_mol) ' GE']);
    %xlabel('Time(s)');
    %ylabel('Intensity (A.U.)');

    
    RE_2=subplot(6,1,3);
    %RE.Position = [.13 .558 .775 .079]; % define your position
       %Simply plotting the  FRET traces now in a subplot below the above plot.
    %A = get(gca,'position');          % gca points at the second one
     %set(gca,'unit','pixels')
    %A(1,4) = A(1,4) * 1.2;              % change the height 
     %A(1,2) = A(1,2) - .5*A(1,4);         % change the vertical position
    %set(gca,'position',A);            % set the values you just changed
    I1=(DonorCorrect_2(i,:));
    I2=((Donor2Correct_2(i,:)));
    I3=(AcceptorCorrect_2(i,:));
    plot(farredlaser_time,I1,'g',farredlaser_time,I2,'r', farredlaser_time,I3,'k');
    ylim([-50 800]);
    %title(['  Molecule ' num2str(n) ' of ' num2str(N_mol) ' GE']);
    %xlabel('Time(s)');
    %ylabel('Intensity (A.U.)');
    
    
    RE_3=subplot(6,1,4);
    %A = get(gca,'position');          % gca points at the second one
    %set(gca,'unit','pixels')
    %A(1,4) = A(1,4) * 1.2;              % change the height 
     %A(1,2) = A(1,2) - .5*A(1,4);         % change the vertical position
    %set(gca,'position',A);            % set the values you just changed
       %Simply plotting the  FRET traces now in a subplot below the above plot.
    I1=(fret01(i,:));
    plot(farredlaser_time,I1,'b');
    ylim([-0.1 1.1]);
    %title(['  Molecule ' num2str(n) ' of ' num2str(N_mol) ' GE']);
    %xlabel('Time(s)');
    %ylabel('Intensity (A.U.)');
    %linkaxes([GE,RE,RE_2],'x')
    
    RE_4=subplot(6,1,5);
     %A = get(gca,'position');          % gca points at the second one
      %set(gca,'unit','pixels')
    %A(1,4) = A(1,4) * 1.2;              % change the height 
     %A(1,2) = A(1,2) - .2*A(1,4);         % change the vertical position
    %set(gca,'position',A);            % set the values you just changed
       %Simply plotting the  FRET traces now in a subplot below the above plot.
    I1=(fret12(i,:));
    plot(farredlaser_time,I1,'b');
    ylim([-0.1 1.1]);
    %title(['  Molecule ' num2str(n) ' of ' num2str(N_mol) ' GE']);
    %xlabel('Time(s)');
    %ylabel('Intensity (A.U.)');
    %linkaxes([GE,RE,RE_2],'x')
    
    RE_5=subplot(6,1,6);
     %A = get(gca,'position');          % gca points at the second one
     %set(gca,'unit','pixels')
    %A(1,4) = A(1,4) * 1.2;              % change the height 
    %A(1,2) = A(1,2) - .2*A(1,4);         % change the vertical position
    %set(gca,'position',A);            % set the values you just changed
       %Simply plotting the  FRET traces now in a subplot below the above plot.
    I1=(fret02(i,:));
    plot(farredlaser_time,I1,'b');
    ylim([-0.1 1.1]);
    %title(['  Molecule ' num2str(n) ' of ' num2str(N_mol) ' GE']);
    %xlabel('Time(s)');
    %ylabel('Intensity (A.U.)');
    %linkaxes([GE,RE,RE_2],'x')
    linkaxes([GE,RE,RE_2,RE_3,RE_4,RE_5],'x')
    
%     Afig=get(hdl,'position');
% Afig(1,4) = Afig(1,4) * 1.5;
% set(hdl,'position',Afig);
%     
    %pause
    
     choice = input('Press 0 to go back one trace\nPress 1 to save current trace\nPress 2 to subtract background (2 clicks)\nPress 3 to cut out\nPress 4 to save the trace as a picture\nPress 5 to skip to data analysis and exit\nPress 6 to collect dwell times\nPress 7 to choose a trace\nPress 8 to count photobleaching steps\nPress 9 to increment counter\nPress Enter to move on\n');

          if choice == 3 
      fprintf('Click twice to select the range to extract\n');
      Done_Cutting_Choice = 0;
      Cut_Counter = 0;
      while Done_Cutting_Choice == 0
           Cut_Counter = Cut_Counter+1;
           [x, y, button] = ginput(2);  % Make Two clicks to specify a region which you want to cut out.
           x(1) = round(x(1)/timeunit); % Starting point of the cut region
           x(2) = round(x(2)/timeunit); % End point of the cut region
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
           
           fid2=num2str(fid2);
           i=num2str(i);
           
           fname1=[pth '\' newfolder '\hel' fid2 ' tr' i '-' num2str(Cut_Counter) '.dat'];
           Output_to_be_saved_inaFile=[fret12(i,x(1):x(2))']
           %Output_to_be_saved_inaFile=[greenlaser_time(x(1):x(2))' DonorCorrect_0(i,x(1):x(2))' Donor2Correct_0(i,x(1):x(2))' AcceptorCorrect_0(i,x(1):x(2))' redlaser_time(x(1):x(2))' DonorCorrect_1(i,x(1):x(2))' Donor2Correct_1(i,x(1):x(2))' AcceptorCorrect_1(i,x(1):x(2))' farredlaser_time(x(1):x(2))' DonorCorrect_2(i,x(1):x(2))' Donor2Correct_2(i,x(1):x(2))' AcceptorCorrect_2(i,x(1):x(2))' fret12(i,x(1):x(2))'];
           save(fname1,'Output_to_be_saved_inaFile','-ascii') ;
           % Asking you whether you want to keep cutting the trace or want
           % to move on to the NEXT set of traces.
           Done_Cutting_Choice=input('\nPress 0 to continue cutting the same trace \nPress Enter to move to next trace\n');   
      end
     end
     
     
     if choice == 9
         [time,y,button]=ginput(1);
         if button==1
        counter_1 = counter_1 + 1;
        fprintf('Current counter = %d\n', counter_1);
         end
         if button==3
                counter_transition=counter_transition + 1;
                fprintf('Current counter transition = %d\n', counter_transition);
         end
         
     end
     if choice ==2 
         %Define X1, X2,X3 and Y1

%GreenExcitation
I1=(DonorCorrect_0(i,:));
I2=((Donor2Correct_0(i,:)));
I3=(AcceptorCorrect_0(i,:));

Ymatrix1=[I1; I2; I3];
X1=greenlaser_time;

%RedExcitation
 I1=(DonorCorrect_1(i,:));
 I2=((Donor2Correct_1(i,:)));
 I3=(AcceptorCorrect_1(i,:));
 
 Ymatrix2=[I1; I2; I3];
 X2=redlaser_time;
 
 %FarRedExcitation
  I1=(DonorCorrect_2(i,:));
  I2=((Donor2Correct_2(i,:)));
  I3=(AcceptorCorrect_2(i,:));
  
  Ymatrix3=[I1; I2; I3];
  X3=farredlaser_time;
  
  %Cy5Cy7FRET
  I1=(fret12(i,:));
  Ymatrix4=[I1];
  X4=X2;
   
     
figurescriptfornucleosomeunwraping(X1, Ymatrix1, X2, Ymatrix2, X3, Ymatrix3, X4, Ymatrix4)
  
     end
    
end
%     subplot(2,1,2);
%     fretE=Cy5(i,:)./(Cy5(i,:)+Cy3(i,:));
%     for m=1:len
%         if Cy5(i,m)+Cy3(i,m)==0
%             fretE(m)=-0.5;
%         end
%     end % This is to avoid undefined fretE interfering with future analysis
%     plot(time,fretE,'b');
%     axis tight;
%     temp=axis;
%     temp(3)=-0.1;
%     temp(4)=1.1;
%     axis(temp);
%     grid on;
%     zoom on;



% close all;
% fclose('all');






