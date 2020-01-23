% pre-process dwelltime.dat files before using Curve Fitting Toolbox to
% fit the dwelltime data from Trace_ViewerO_smarcal.m
% 
% use together with dwelltime_fit2plot.m
%
% by X. Feng xfeng17@jhu.edu

% To fit the dwelltime histogram,
% 0. input your filename below, then run the script
% 1. go to APPS > Curve Fitting
% 2. select 'centers' as x and 'hist_values' as y
% 3. fit the data with a function of your choice
% 4. export the fitting results: Fit > Save to workspace > OK
% 5. Use dwelltime_fit2plot.m to plot the fit and the histogram

% USER INPUT
filename = 'dwelltime1_0806-0820.dat';

% End of user input

dt = readtable(filename);
dwelltimes = dt.Var1;

[hist_values, edges] = histcounts(dwelltimes);

n_bins = length(edges) - 1;
bin_width = edges(2) - edges(1);

centers = edges(1:n_bins) + bin_width/2;


