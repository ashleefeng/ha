% plot percent observations
listing = dir(pwd);
fig = figure;
files = {listing.name};
i = 1;
n_good = 0;
n_bad = 0;
n_hel1_exch = 0; %127 for + ATP, 2 for -ATP, 1 for ATPgS

while i <= length(files)
    fname = files{i};
    if strcmp(fname,'.') || strcmp(fname,'..')
        i = i + 1;
        continue
    end
    
    [filepath,name,ext] = fileparts(fname);
   	if strcmp(ext, '.dat')
        fprintf('Trace %d / %d %s\n', i - 2, length(files) - 2, fname);
        opts = detectImportOptions(fname);
        opts.VariableNames = {'trcID', 'is_junk'};
        junk_counts = readtable(fname, opts);
        n_bad = n_bad + sum(junk_counts.is_junk);
        n_good = n_good + height(junk_counts) - sum(junk_counts.is_junk);
        disp(1 - sum(junk_counts.is_junk) / height(junk_counts));
    end
    i = i + 1;
end

fprintf('Fraction of traces that show exchange in nuc only: %.1f%% (%d / %d)\n', n_hel1_exch/n_good*100, n_hel1_exch, n_good);