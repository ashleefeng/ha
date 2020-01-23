% plot the fit from Curve Fittig Toolbox 
%
% See dwelltime_dat2fit.m for usage
%
% by X. Feng xfeng17@jhu.edu

model = fittedmodel;
goodness = goodness;
f = figure();
set(gcf,'Position',[500 500 600 400])
bar(centers, hist_values, 1.0, 'FaceColor', [0.3010 0.7450 0.9330]);

hold on

xmin = min(centers);
xmax = max(centers);
step = (xmax-xmin)/100;
x = (xmin):step:(xmax+step);
fit = plot(x, model(x), 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Count');

model
goodness
