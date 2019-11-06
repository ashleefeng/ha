% view 3 color traces exported as .dat files
% X. Feng Nov 6, 2019 xfeng17@jhu.edu

listing = dir(pwd);
figure;
for file = {listing.name}
    fname = file{1};
    if strcmp(fname,'.') || strcmp(fname,'..')
        continue
    end
    [filepath,name,ext] = fileparts(fname);
   	if strcmp(ext, '.dat')
        ginput;
        trc_data = three_color_trace_dat_reader(fname);
        three_color_trace_plotter(trc_data, name, 1, 500);
    end
end