traceID2count = readmatrix('bound_traces.csv');
n = size(traceID2count, 1);
counts = traceID2count(:, 2);
s = sum(counts);
fraction_bound = s/n;
disp(s);
disp(n);
disp(fraction_bound);
