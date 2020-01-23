dt = dt4;

h = histogram(dt, 10);
xlabel('Time(s)');
ylabel('Count');
set(gca, 'FontSize', 15);

mean(h.Values)
h.BinEdges
x = h.BinEdges + (h.BinEdges(2) - h.BinEdges(1))/2;
x = x(1:size(x,2)-1);

y = h.Values;

xfit = 0.9 :0.1:20;
hold on
plot(xfit, fittedmodel(xfit), 'LineStyle', '--', 'LineWidth', 2);
title('\tau_4 = 1.63s (1.37, 2.02)');
% xlim([-10 200]);
% ylim([0 170]);
hold off
diary dt4_fit.log
fittedmodel
goodness
diary off