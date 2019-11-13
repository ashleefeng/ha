% view 3 color traces exported as .dat files
% X. Feng Nov 6, 2019 xfeng17@jhu.edu

listing = dir(pwd);
figure;
files = {listing.name};
i = 3;
while i < length(files)
    fname = files{i};
    if strcmp(fname,'.') || strcmp(fname,'..')
        continue
    end
    
    [filepath,name,ext] = fileparts(fname);
   	if strcmp(ext, '.dat')
        keyanswer = input('Enter=next trace, b=go back: ', 's');
        if strcmp(keyanswer, 'b')
            i = i - 1;
        end
        trc_data = three_color_trace_dat_reader(fname);
        three_color_trace_plotter(trc_data, name, 1, 500);
    end
    i = i + 1;
end