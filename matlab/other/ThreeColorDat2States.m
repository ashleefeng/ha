
% Read .dat file of a selected trace and convert it into a color-coded 
%   line representating state transitions over time.
% Requires utility functions: which_state, bd_filter, and getDT.
% Jan 8, 2023 X. Feng

close all
clear all

% ---------------------------- USER INPUT ----------------------------

trace_prefix = 'hel1_trace_40';
num_std = 5;
save_output = 1; % 1 if saving output; 0 if not saving output
time_unit = 0.05; % unit is seconds

% -------------------------- USER INPUT ENDS --------------------------

% Read input trace data

trace_fname = [trace_prefix '.dat'];
data = readmatrix(trace_fname);

background_fname = [trace_prefix '_background.dat'];
opts = detectImportOptions(background_fname, 'LeadingDelimitersRule', 'ignore');
background = readmatrix(background_fname, opts);

x = data(:, 1) ;
y0 = data(:, 2);
y1 = data(:, 3);
y2 = data(:, 4);
u0 = background(4);
u1 = background(5);
u2 = background(6);

% Plot trace

figure
plot(x, y0, 'g');
hold on
plot(x, y1, 'r');
plot(x, y2, 'b');
%plot(x, y0 + y1 + y2, 'k');
hold off
xlabel("Time (s)");
ylabel("Fluor. Intensity");
set(gcf, 'Position', [500 600 800 100]);

if save_output
    saveas(gcf, trace_prefix, 'epsc');
end

% Assign a state - 0, 1, 2, or 3 - to each data point

y_state = zeros(1, length(x));
    
ptr = 0;
for t = x'
    ptr = ptr + 1; 
    state = which_state(y0(ptr), y1(ptr), y2(ptr), u0, u1, u2, num_std);
    y_state(ptr) = state;        
end

% Create a matlab colormap with 3 or 4 colors, depending on # states in a trace

if ismember(0, y_state)
    cm = [0 1 0; 1 0 0; 0 0 1; 1 1 1];
else
    cm = [1 0 0; 0 0 1; 1 1 1];
end

% Make a heatmap trace of assigned states

figure
colormap(cm);
imagesc(y_state);
set(gcf, 'Position', [500 450 800 50]);

if save_output
    heatmap_fname = [trace_prefix '_states'];
    saveas(gcf, heatmap_fname, 'epsc');
end

% Filter out false positive states outside of binding event

filtered_state_list = bd_filter(y_state);

if ismember(0, filtered_state_list)
    cm = [0 1 0; 1 0 0; 0 0 1; 1 1 1];
else
    cm = [1 0 0; 0 0 1; 1 1 1];
end

figure
colormap(cm);
imagesc(filtered_state_list);
set(gcf, 'Position', [500 300 800 50]);

if save_output
    heatmap_fname = [trace_prefix '_states_filtered'];
    saveas(gcf, heatmap_fname, 'epsc');
end


% Extract dwell times from state list

dt0_list = getDT(filtered_state_list, 0, time_unit);
dt1_list = getDT(filtered_state_list, 1, time_unit);
dt2_list = getDT(filtered_state_list, 2, time_unit);

% Save dwell time data

if save_output
    dt0_fname = [trace_prefix '_cy3onlyDT.dat'];
    dt1_fname = [trace_prefix '_cy5fretDT.dat'];
    dt2_fname = [trace_prefix '_cy7fretDT.dat'];

    csvwrite(dt0_fname, dt0_list);
    csvwrite(dt1_fname, dt1_list);
    csvwrite(dt2_fname, dt2_list);
end
