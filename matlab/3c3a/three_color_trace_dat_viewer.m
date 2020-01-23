% view 3 color traces exported as .dat files
% X. Feng Nov 6, 2019 xfeng17@jhu.edu

listing = dir(pwd);
fig = figure;
files = {listing.name};
i = 3;
fret1start = [];
fret1end = [];
fret1val = [];
fret1id = [];
fret2start = [];
fret2end = [];
fret2val = [];
fret2id = [];
dt3id = [];
dt3start = [];
dt3end = [];
dt3 = [];
fret301 = [];
fret302 = [];
file_prefix = '../analysis/';

while i <= length(files)
    fname = files{i};
    if strcmp(fname,'.') || strcmp(fname,'..')
        continue
    end
    
    [filepath,name,ext] = fileparts(fname);
   	if strcmp(ext, '.dat')
        fprintf('Trace %d / %d %s\n', i - 2, length(files) - 2, fname);
        trc_data = three_color_trace_dat_reader(fname);
        three_color_trace_plotter(trc_data, name, 1, 500);
        keyanswer = input('Enter=next trace, b=go back, f=collect avg fret, g=go to, s=save as png t=terminate: ', 's');
        if strcmp(keyanswer, 'b')
            i = i - 2;
        elseif strcmp(keyanswer, 'f')
            [time,~,button]=ginput;
            [fret1id, fret1start, fret1end, fret1val, ...
                fret2id, fret2start, fret2end, fret2val, ...
                dt3id, dt3start, dt3end, dt3, fret301, fret302] = ...
                dwelltime_collector(trc_data, i-2, time, button, ...
                fret1id, fret1start, fret1end, fret1val, ...
                fret2id, fret2start, fret2end, fret2val, ...
                dt3id, dt3start, dt3end, dt3, fret301, fret302);
        elseif strcmp(keyanswer, 'g')
            i = input('Trace to go to: ') + 1;
        elseif strcmp(keyanswer, 's')
            saveas(fig, ['../plots/' name '.png']);
        elseif strcmp(keyanswer, 't')
            break
        end
        
    end
    i = i + 1;
end

if ~isempty(fret1val)
    fret1 = [fret1id; fret1start; fret1end; fret1val]';
    fname = [file_prefix 'dt1_productive_unwrapping_redo.dat'];
    save(fname,'fret1','-ascii','-append');
end

% if ~isempty(fret2val)
%     fret2 = [fret2id; fret2start; fret2end; fret2val]';
%     fname = [file_prefix '_time2second_exchange.dat'];
%     save(fname,'fret2','-ascii','-append');
% end

if ~isempty(dt3)
    dt3save = [dt3id; dt3start; dt3end; dt3; fret301; fret302]';
    fname = [file_prefix 'time2second_exchange.dat'];
    save(fname,'dt3save','-ascii','-append');
end