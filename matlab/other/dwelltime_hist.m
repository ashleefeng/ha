dt = [dt3a;dt3b];

h = histogram(dt);
xlabel('Time(s)');
ylabel('Count');
set(gca, 'FontSize', 20);

mean(h.Values)
h.BinEdges
x = h.BinEdges + (h.BinEdges(2) - h.BinEdges(1))/2;
x = x(1:size(x,2)-1);

y = h.Values;

xfit = 10:1:260;
hold on
plot(xfit, fittedmodel(xfit), 'LineStyle', '--', 'LineWidth', 2);
title('\tau_3 = 48.2 s (38.7, 63.9)');