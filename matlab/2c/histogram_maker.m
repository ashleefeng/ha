% plot intensity and fret histograms from .traces 
% by X. A. Feng
% Feb 23, 2019

% USER INPUT
FRAME_RATE = 0.1; % s/frame
LEAKAGE = 0.12;
GAMMA = 1;
RANGE = [[12 16]; [2, 6]]; % Put fret frames first!

% filters
DON_ACC_SUM_RANGE = [0,2000]; % filter for donor + acceptor intensity
ACC_RANGE = [300, 1500]; % filter for acceptor intensity
FRET_RANGE = [-1.2, 1.2]; % filter for FRET efficiency
% End of user input

directory = input('Directory: ','s');

if isempty(directory)
    directory = pwd;
end

cd(directory);

trace_files = dir('*.traces');

n_traces_files = size(trace_files, 1);

n_molecs = 0;

for i = 1:n_traces_files
    trace_file = fopen(trace_files(i).name);
    n_frames=fread(trace_file, 1, 'int32');% This fetches the total duration of the imaging that was carried out for the
    
    
    num_traces=fread(trace_file, 1, 'int16');
    n_molecs = n_molecs + num_traces / 2;
    
    %disp(strcat('The number of traces in this file is: ', string(num_traces / 2)));
    fclose(trace_file);

end

don_mean_all = zeros(n_molecs, 2);
acc_mean_all = zeros(n_molecs, 2);
fret_all = zeros(n_molecs, 2);

molec_id = 1;

for i = 1:n_traces_files
    trace_file = fopen(trace_files(i).name);
    n_frames=fread(trace_file, 1, 'int32');
    
    num_traces=fread(trace_file, 1, 'int16');
    num_molecules = num_traces / 2;
    
    raw_data=fread(trace_file, num_traces*n_frames, 'int16');
    fclose(trace_file);
    
    indices = (1: num_traces * n_frames);
    data = zeros(num_traces, n_frames);
    don = zeros(num_traces / 2, n_frames);
    acc = zeros(num_traces / 2, n_frames);
    data(indices) = raw_data(indices);
    
    for j = 1: (num_traces/2)
        don(j, :)= data(j * 2  - 1, :);
        acc(j, :)= GAMMA .* data(j * 2, :);
    
        for k = 1: 2
            foi = RANGE(k, 1): RANGE(k, 2);
            don_mean_all(molec_id, k) = mean(don(j, foi), 2);
            acc_mean_all(molec_id, k) = mean(acc(j, foi), 2);
            fret_all(molec_id, k) = mean((acc(j, foi)...
                                          - LEAKAGE * don(j, foi))...
                                      ./ (acc(j, foi)...
                                          - LEAKAGE * don(j, foi)...
                                          + (don(j, foi))), 2);
        end
        molec_id = molec_id + 1;
    end
end

n_filtered = 0;

for i = 1: n_molecs
    da_sum = don_mean_all(i, 1) + acc_mean_all(i, 1);
    acc_mean = acc_mean_all(i, 2);
    fret = fret_all(i, 1);
    
    if da_sum < DON_ACC_SUM_RANGE(1) || da_sum > DON_ACC_SUM_RANGE(2) ||...
        acc_mean < ACC_RANGE(1) || acc_mean > ACC_RANGE(2) || ...
        fret < FRET_RANGE(1) || fret > FRET_RANGE(2)
        
        fret_all(i) = nan;
        don_mean_all(i, :) = [nan nan];
        acc_mean_all(i, :) = [nan nan];
    else
        n_filtered = n_filtered + 1;
    end
    
end

fprintf("The plots include %d out of %d (%.2f %%) total traces. \n", n_filtered, n_molecs, n_filtered/n_molecs*100);

fig = figure('Position',[100, 100, 1400, 600]);
for k = 1:2
    subplot(2, 4, (k-1)*4 + 1);
    histogram(don_mean_all(:, k) + acc_mean_all(:, k), 100, 'FaceColor', 'k', 'EdgeColor', [0 0 0]);
    title(strcat("Frames ", string(RANGE(k, 1)), " to ", string(RANGE(k, 2)), " donor + acceptor"));
    xlabel('Intensity');
    ylabel('Count');
    
    subplot(2, 4, (k-1)*4 + 2);
    histogram(don_mean_all(:, k), 100, 'FaceColor', [0 0.5 0] , 'EdgeColor',  [0 0.5 0]);
    title(strcat("Frames ", string(RANGE(k, 1)), " to ", string(RANGE(k, 2)), " donor"));
    xlabel('Intensity');
    ylabel('Count');
    
    subplot(2, 4, (k-1)*4 + 3);
    histogram(acc_mean_all(:, k), 100, 'FaceColor', 'r', 'EdgeColor', 'r');
    title(strcat("Frames ", string(RANGE(k, 1)), " to ", string(RANGE(k, 2)), " acceptor"));
    xlabel('Intensity');
    ylabel('Count');
    
    subplot(2, 4, (k-1)*4 + 4);
    histogram(fret_all(:, k), 'BinWidth', 0.02, 'FaceColor', 'b', 'EdgeColor', 'b');
    title(strcat("Frames ", string(RANGE(k, 1)), " to ", string(RANGE(k, 2)), " FRET"));
    xlim([-0.1, 1.1]);
    xlabel('FRET efficiency');
    ylabel('Count');
end

savefig('histograms_190228b.fig');
saveas(fig, 'histograms_190228b.png');

% Save FRET values
fret_list = fret_all(:, 1);
to_save = fret_list(~isnan(fret_list));
csvwrite('fret_values.csv', to_save);