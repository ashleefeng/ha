%monte carlo simulation for nucleosome unwrapping simultaneously occuring with
% h2a eviction'three color data

%Matt Poyton September 20, 2019

%first we need an equation for the photobleaching lifetime of Cy3, the
%frequency of low-fret, dna unwrapping events and the photobleaching rates
%of cy5 and cy7. For the sake of simplicity, let's say that Cy3 has a
%lifetime of 60 seconds, Cy5 has a lifetime of 45 seconds and Cy7 has a
%lifetime of 30 seconds. The nucleosome will be unwrapped for
% 5% of the total time of the trace, which will be
% set to 300 seconds. Now we need the cumulative
%distributions for all of these lifetimes.

Cy3Lifetime = 60;
Cy5Lifetime = 45;
Cy7Lifetime = 30;
%UnwrappingEquilibrium = 0.05;
Cy3Increase = 0; %If a randomly generated number is greater than 0 than Cy3 will increase
number_unwrapping_events = 20; %nuber of unwrapping events per trace

n_trials = 1E6; %number of random trials
ObservationTime=300; %

PredictedCy3BleachingTime=zeros(n_trials,1);
PredictedCy5BleachingTime=zeros(n_trials,1);
PredictedCy7BleachingTime=zeros(n_trials,1);
Unwrapped=zeros(n_trials,10);
DiditColocalize=zeros(n_trials,number_unwrapping_events);
Cy3IncreaseProb=NaN(n_trials,1);

%This section generates 3 random lifetimes for Cy3, Cy5 and Cy7 using
%the initialized variables above

for i=1:n_trials
    y=rand(1);
PredictedCy3BleachingTime(i)=log(y)*Cy3Lifetime*-1;
end

for i=1:n_trials
    y=rand(1);
PredictedCy5BleachingTime(i)=log(y)*Cy5Lifetime*-1;
end

for i=1:n_trials
    y=rand(1);
PredictedCy7BleachingTime(i)=log(y)*Cy7Lifetime*-1;
end

%Ok, we generate the time at which the 'number_unwrapping_events' each
%occur

for i=1:n_trials
    for j=1:number_unwrapping_events
    %y=rand(1);
%if y < UnwrappingEquilibrium;
    Unwrapped(i,j)=rand(1)*ObservationTime;
%else
   % Unwrapped(i,n)=NaN;
end
    end
%end

%We generate the one time in the trace at which Cy3 increases in intensity

for i=1:n_trials
    y=rand(1);
if y > Cy3Increase;
    Cy3IncreaseProb(i)=rand(1)*ObservationTime;
else
    Cy3IncreaseProb(i)=NaN;
end
end

%If the unwrapping event occurs within +/-0.6 seconds of the unwrapping
%event and before any of the dyes photobleach that is a positive event.


for i=1:n_trials
    for j=1:number_unwrapping_events
if Cy3IncreaseProb(i)> Unwrapped(i,j) -.6 & Cy3IncreaseProb(i) < Unwrapped(i,j) +.6  & Cy3IncreaseProb(i)< PredictedCy3BleachingTime(i) & Cy3IncreaseProb(i)< PredictedCy5BleachingTime(i) & Cy3IncreaseProb(i)<PredictedCy7BleachingTime(i);
    DiditColocalize(i,j)=1;
    else
    DiditColocalize(i,j)=0;
end
    end
end

%The fraction of traces that show unwrapping events

DiditColocalize2=sum(DiditColocalize,2);

counter=sum(DiditColocalize2);

finalfraction = counter/n_trials;

display(finalfraction)
