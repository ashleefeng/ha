clear
hel = csvread('hel1hel2_photobleaching.csv');
n = size(hel, 1);
for i = 1:n
    if hel(i) == -1
        hel(i) = nan;
    end
end

[hist_values, edges] = histcounts(hel);

n_bins = length(edges) - 1;
bin_width = edges(2) - edges(1);

centers = edges(1:n_bins) + bin_width/2;