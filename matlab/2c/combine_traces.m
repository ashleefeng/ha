%%% conversion of selected traces into .dat for ebFRET
%%% put all selected traces into a folder and run
%%% David Gnutt 11/12/2016


Directory_of_Traces=input('Directory:: ');
    if isempty(Directory_of_Traces)
       Directory_of_Traces=pwd;
    end
   
cd(Directory_of_Traces);
   
list_of_files=dir([fullfile(Directory_of_Traces) '\\*.dat']); 
len=length(list_of_files);

%load traces into matlab
for i = 1:len
    Traces{i}=importdata(list_of_files(i).name);
    len_traces{i}=length(Traces{i});
end
   
b=cell2mat(len_traces);
%%% create a stacked set of traces ---- each trace gets a unique
%%% identifier, i.e. 1, 2, 3 ... N in column 1
%%% column 2: D(trace1,time1),...D(trace1,timepoint N),...
%%%           D(trace N,time N)
%%% column 3: A(trace1, time1),...A(trace1, time N),...
%%%           A(trace N,time N)
for i= 1:len 
      
       if i == 1
           stacked_traces(1:len_traces{i},2:3)=Traces{i}(:,2:3);
           stacked_traces(1:len_traces{i},1)=i;
       else
           stacked_traces(1+sum(b(1:i-1)):sum(b(1:i)),2:3)=Traces{i}(:,2:3);
           stacked_traces(1+sum(b(1:i-1)):sum(b(1:i)),1)=i;
       end
end
%%% saves file as ASCII .dat files with 3 columns - this file can be loaded into "ebFRET" and analyzed.   
save('stacked_traces.dat', 'stacked_traces', '-ascii');