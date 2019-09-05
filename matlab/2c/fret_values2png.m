fret_values = readtable('fret_values.csv', 'ReadVariableNames', false);
fig = histogram(fret_values.Var1);
set(gcf,'Position',[100 100 600 400]);
set(gca, 'FontSize', 20)
xlim([-0.2 1.2]);
xlabel('FRET Efficiency');
ylabel('Count');
saveas(fig, 'matlab_histogram.png');