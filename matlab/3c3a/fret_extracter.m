% view 3 color traces exported as .dat files
% X. Feng Nov 6, 2019 xfeng17@jhu.edu

listing = dir(pwd);
files = {listing.name};
i = 1;
fret1start = [];
fret1end = [];
fret1val = [];
fret1id = [];
fret2start = [];
fret2end = [];
fret2val = [];
fret2id = [];
file_prefix = 'cy3cy5';

while i <= length(fret2all.e00)
    traceID = fret2all.e00(i);
    fname = files{traceID + 2};
    
    if strcmp(fname,'.') || strcmp(fname,'..')
        continue
    end
    
    [filepath,name,ext] = fileparts(fname);
    
    if strcmp(ext, '.dat')
        fprintf('Trace %d / %d\n', traceID, length(files) - 2);
        trc_data = three_color_trace_dat_reader(fname);
        time = [fret2all.e01(i) fret2all.e1(i)];
        button = [2 2];
        [fret1id, fret1start, fret1end, fret1val, ...
            fret2id, fret2start, fret2end, fret2val] = ...
            dwelltime_collector(trc_data, traceID, time, button, ...
            fret1id, fret1start, fret1end, fret1val, ...
            fret2id, fret2start, fret2end, fret2val);
    end
    
    i = i + 1;
end

% if ~isempty(fret1val)
%     fret1 = [fret1id; fret1start; fret1end; fret1val]';
%     fname = [file_prefix '_fret1_unwrapping.dat'];
%     save(fname,'fret1','-ascii','-append');
% end

if ~isempty(fret2val)
    fret2 = [fret2id; fret2start; fret2end; fret2val]';
    fname = [file_prefix '_fret2_static.dat'];
    save(fname,'fret2','-ascii','-append');
end