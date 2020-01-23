values = dt4;

[hist_values, edges] = histcounts(values, 10);

n_bins = length(edges) - 1;
bin_width = edges(2) - edges(1);

centers = edges(1:n_bins) + bin_width/2;
