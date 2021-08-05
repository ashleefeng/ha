function fret_values = fret_values_extractor(trc, time, timeunit, button)
    time1 = time(button == 1);
    t1 = ceil(time1(1) / (3 * timeunit));
    t2 = ceil(time1(2) / (3 * timeunit));
    fret_values = trc.fret01(t1:t2);
end