DIR = dir;
set(0,'defaultfigurewindowstyle','docked');
fig = figure;
j = 3;
while true
    FILENAME = DIR(j).name;
    hel = readtable(FILENAME);
    time = table2array(hel(:, 1));
    don = table2array(hel(:, 2));
    acc = table2array(hel(:, 3));
    sum = don + acc + ones(1000, 1) * 500;
    fret= acc./(acc + don);
    for i = 1:1000
        if ((fret(i, 1) < -0.2) || (fret(i, 1) > 1.2))
            fret(i, 1) = NaN;
        end
    end
    subplot(2, 1, 1);
    
    plot(time, don, 'g');
    hold on
    plot(time, acc, 'r');
    plot(time, sum, 'b');
    xlabel('Time(s)');
    ylabel('Intensity (a. u.)');
    title(FILENAME);
    hold off
    subplot(2, 1, 2);
    plot(time, fret, 'b');
    xlabel('Time(s)');
    ylabel('FRET efficiency');
    figure(fig);
    in = input('Enter q to exit, b to go back\n', 's');
    if in == 'q'
        break;
    end
    if in == 'b'
        j = j - 2;
    end
    j = j + 1;
    if j > size(DIR, 1)
        j = j - 1;
        disp('This is the last trace!');
    end
    if j < 3
        j = 3;
    end
end

close();