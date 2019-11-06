dt = g1;

h = histogram(dt, 'BinWidth', 5);
xlabel('Time(s)');
ylabel('Count');
set(gca, 'FontSize', 20);

mean(h.Values)
h.BinEdges
x = h.BinEdges + (h.BinEdges(2) - h.BinEdges(1))/2;
x = x(1:size(x,2)-1);

y = h.Values;

xfit = 2.5:1:200;
hold on
plot(xfit, fittedmodel(xfit), 'LineStyle', '--', 'LineWidth', 2);
title('\tau_2 = 8.85s (8.18, 9.65)');
xlim([-10 200]);
ylim([0 170]);
hold off
diary ATPgS_fit.log
fittedmodel
goodness
diary off