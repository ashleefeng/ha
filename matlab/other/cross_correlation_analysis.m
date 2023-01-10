% Calculate pearson correlation coeff for time-series pairs
close all

folder = '/Users/ashleefeng/OneDrive - Johns Hopkins/Lab-Data Backup/NURF_and_GAGA/Analysis/462_DNA_anticorrelation/raw_data';
file_pattern = fullfile(folder, '**/*.dat');
files = dir(file_pattern);
n_files = length(files);
pearson_r_list = zeros(n_files, 1);
lower_r_list = zeros(n_files, 1); % 95% confidence interval lower bound
upper_r_list = zeros(n_files, 1);% 95% confidence interval upper bound

figure

for i = 1 : n_files
    filename = files(i).name;
    fprintf(1, 'Now reading file #%d: %s\n', i, filename);
    full_filename = fullfile(files(i).folder, filename);
    data = readtable(full_filename);
    data = rmmissing(data);
    data.Properties.VariableNames = {'time' 'fret35' 'fret37'};
    y1 = data.fret35;
    y2 = data.fret37;
    if (i == 15) || (i == 23) || (i == 42) || (i == 24) || (i == 38)
        
        scatter(y1, y2, 'o', 'LineWidth', 2)
        
        hold on
    end
    [R, P, RL, RU] = corrcoef(y1, y2);
    pearson_r_list(i) = R(1, 2);
    lower_r_list(i) = RL(1, 2);
    upper_r_list(i) = RU(1, 2);
end

hold off
xlabel('Cy3-Cy5 FRET');
ylabel('Cy3-Cy7 FRET');
legend('Trace 1', 'Trace 2', 'Trace 3', 'Trace 4', 'Trace 5');
set(gca, 'FontSize', 25);

figure
x = 1:n_files;
[out, idx] = sort(pearson_r_list);

xconf = [x x(end:-1:1)];
upper_sorted = upper_r_list(idx);
lower_sorted = lower_r_list(idx);
yconf = [upper_sorted; lower_sorted(end:-1:1)]';
p = fill(xconf, yconf, 'red');
p.EdgeColor = 'none';
p.FaceColor = [1 0.8 0.8];
hold on
plot(x, pearson_r_list(idx), 'or', 'LineWidth', 2);
plot(x, zeros(n_files, 1), ':k', 'LineWidth', 3);
set(gca, 'FontSize', 25);
ylim([-1 1]);
xlim([0 70]);
xlabel('Trace # Sorted by r');
ylabel('Pearson r');
legend('95% confidence interval', 'Pearson r');
hold off