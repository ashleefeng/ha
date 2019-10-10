values = csvread('Cy3Cy5_fret_values_ver2.csv');

[hist_values, edges] = histcounts(values, 'BinWidth', 0.025);

n_bins = length(edges) - 1;
bin_width = edges(2) - edges(1);

centers = edges(1:n_bins) + bin_width/2;
