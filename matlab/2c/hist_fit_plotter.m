FILENAME = 'fret.dat';
TITLE = 'Nucleosome only';
SAVEAS = 'fret_hist_matlab_180117';

frethist = readtable(FILENAME, 'Delimiter', '\t');

x = table2array(frethist(:, 1));
x = x(1:size(x, 1));
y = table2array(frethist(:, 3));
y = y(1:size(y, 1));

f = bar(x, y);
f.FaceColor = 'b';
f.FaceAlpha = 0.1;
f.EdgeColor = 'b';
f.LineWidth = 2;
f.BarWidth = 1;
% title(TITLE);
% hold on
% 
% fitx = table2array(fit(:, 1));
% fit1 = table2array(fit(:, 3));
% fit2 = table2array(fit(:, 4));
% fit3 = table2array(fit(:, 5));
% 
% plot(fitx, fit1, '--', 'LineWidth', 3, 'color', 'm');
% plot(fitx, fit2, '--', 'LineWidth', 3, 'color', 'g');
% plot(fitx, fit3, '--', 'LineWidth', 3, 'color', 'r');

xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 50);
xlabel('FRET efficiency', 'FontSize', 50);
ylabel('Frequency', 'FontSize', 50);
xlim([0, 1]);
ylim([0, 0.05]);
saveas(gcf, strcat(SAVEAS, '.png'));
% saveas(gcf, strcat(SAVEAS, '.fig'));