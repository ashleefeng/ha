function [fret1id, fret1start, fret1end, fret1val, fret2id, fret2start, fret2end, fret2val, dt3id, dt3start, dt3end, dt3, fret301, fret302] = dwelltime_collector(trc, id, time, button, fret1id, fret1start, fret1end, fret1val, fret2id, fret2start, fret2end, fret2val, dt3id, dt3start, dt3end, dt3, fret301, fret302)
    
    timeunit = 0.2;

    time1 = time(button == 1);
    
    for c = 1:2:(sum(button == 1) - 1)
        
        t1 = ceil(time1(c) / (3 * timeunit));
        t2 = ceil(time1(c+1) / (3 * timeunit));
        
        start_pt = time1(c);
        end_pt = time1(c+1);
        
        fret1id(end + 1) = id;
        fret1start(end + 1) = start_pt;
        fret1end(end + 1) = end_pt;
        fret1val(end + 1) = mean(trc.fret12(t1:t2));
        
    end

    time2 = time(button == 2);
    
    for c = 1:2:(sum(button == 2) - 1)
        
        t1 = ceil(time2(c) / (3 * timeunit));
        t2 = ceil(time2(c+1) / (3 * timeunit));
        
        start_pt = time2(c);
        end_pt = time2(c+1);
        
        fret2id(end + 1) = id;
        fret2start(end + 1) = start_pt;
        fret2end(end + 1) = end_pt;
        fret2val(end + 1) = mean(trc.fret12(t1:t2));
        
    end
    
    time3 = time(button == 3);
    
    for c = 1:2:(sum(button == 3) - 1)
        
        t1 = ceil(time3(c) / (3 * timeunit));
        t2 = ceil(time3(c+1) / (3 * timeunit));
        
        start_pt = time3(c);
        end_pt = time3(c+1);
        
        dt3id(end + 1) = id;
        dt3start(end + 1) = start_pt;
        dt3end(end + 1) = end_pt;
        dt3(end + 1) = abs(end_pt - start_pt);
        fret301(end + 1) = mean(trc.fret01(t1:t2));
        fret302(end + 1) = mean(trc.fret02(t1:t2));
        
    end
    
end