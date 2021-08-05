files = dir('*.dat');
n_files = size(files, 1);
total_length = 0;

for i = 1: n_files
    filename = files(i).name;
    subtrace = readtable(filename);
    subtrace_length = size(subtrace, 1);
    total_length = total_length + subtrace_length;
    %break
end

avg_fret_list = zeros(n_files, 1);
all_fret_list = zeros(total_length, 3);
ptr = 1;

for i = 1: n_files
    filename = files(i).name;
    subtrace = readtable(filename);
    subtrace_length = size(subtrace, 1);
    subtrace_fret = subtrace.Var3./(subtrace.Var2 + subtrace.Var3);
    all_fret_list(ptr:(ptr + subtrace_length - 1), 1) = i;
    all_fret_list(ptr:(ptr + subtrace_length - 1), 2) = subtrace_fret;
    all_fret_list(ptr:(ptr + subtrace_length - 1), 3) = 1/subtrace_length;
    ptr = ptr + subtrace_length;
    
    avg_fret = mean(subtrace_fret);
    if (avg_fret < -0.2) || (avg_fret > 1.1)
        avg_fret_list(i) = NaN;
    else
        avg_fret_list(i) = avg_fret;
    end
end

for j = 1 : total_length
    fret_val = all_fret_list(j, 2);
    if (fret_val < -0.1) || (fret_val > 1.2)
        all_fret_list(j, 2) = NaN;
    end
end

[hist, edges] = histwc(all_fret_list(:, 2), all_fret_list(:, 3), 80);

figure
bar(edges, hist);
xlim([-0.2 1.1]);
xlabel('FRET Efficiency');
ylabel('Weighted Count');
set(gca, 'fontsize', 25);
