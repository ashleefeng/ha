function three_color_trace_plotter(trc, trc_name, firstpoint, lastpoint)
    
    subplot('position',[0.1 0.84 0.8 0.12]); 
    
    plot(trc.time, trc.ex0em0, 'g', trc.time, trc.ex0em1, 'r', trc.time, trc.ex0em2, 'm');
    hold on
    plot(trc.time, trc.ex0em0 + trc.ex0em1 + trc.ex0em2 + 200, 'k');
    hold off
    title(['Green Laser Molecule ' trc_name], 'Interpreter', 'none');
    %title('Green laser on');
    %ylabel('Intensity');
    zoom on;
    set(gca, 'FontSize', 12);
    %ylim([-100, 500]);
    set(gca,'XTickLabel',[]);
    
    % red laser excitation corrected trace
    subplot('position',[0.1 0.68 0.8 0.12]); 
    %plot(trc.time, trc.ex1em0, 'g', trc.time, trc.ex1em1, 'r', trc.time, trc.ex1em2, 'm');
    plot(trc.time, trc.ex1em1, 'r', trc.time, trc.ex1em2, 'm');
    title(['Red Laser Molecule ' trc_name], 'Interpreter', 'none');
    %title('Red laser on');
    zoom on;
    %ylabel('Intensity');
    set(gca, 'FontSize', 12);
    %ylim([-100, 600]);
    set(gca,'XTickLabel',[]);
    
    % 750 laser excitation corrected trace
    subplot('position',[0.1 0.52 0.8 0.12]);
    %plot(trc.time, trc.ex2em0, 'g', trc.time, trc.ex2em1, 'r', trc.time, trc.ex2em2, 'm');
    plot(trc.time, trc.ex2em2, 'm');
    title(['750 Laser Molecule ' trc_name], 'Interpreter', 'none');
    %title('750 laser on');
    zoom on;
    %ylabel('Intensity');
    set(gca, 'FontSize', 12);
    %ylim([-50, 250]);
    %xlabel('Time (s)');
    
    subplot('position',[0.1 0.36 0.8 0.1]);
    %	FretEc=(1./(1+gamma*(donorcorrect(i,:)./acceptorcorrect(i,:))));
    plot(trc.time(firstpoint:lastpoint), trc.fret01(firstpoint:lastpoint));
    ylim([-0.1 1.1]);
    title('Cy3 Cy5 FRET');
    ylabel('FRET');
    set(gca, 'FontSize', 12);
    
    subplot('position',[0.93 0.36 0.03 0.1]);
    x = -0.1:0.05:1.1;
    [hX,hN]=hist(trc.fret01(firstpoint:lastpoint),x);
    barh(hN,hX,'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    xlabel('Count');
    grid on;
    axis on;
    zoom on;
    
    subplot('position',[0.1 0.20 0.8 0.1]);
    plot(trc.time(firstpoint:lastpoint), trc.fret02(firstpoint:lastpoint));
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    zoom on;
    title('Cy3 Cy7 FRET');
    ylabel('FRET');
    set(gca, 'FontSize', 12);
    
    subplot('position',[0.93 0.20 0.03 0.1]);
    x = -0.1:0.05:1.1;
    [hX,hN]=hist(trc.fret02(firstpoint:lastpoint),x);
    barh(hN,hX,'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    xlabel('Count');
    grid on;
    axis on;
    zoom on;
    
    subplot('position', [0.1 0.04 0.8 0.1]);
    plot(trc.time(firstpoint:lastpoint), trc.fret12(firstpoint:lastpoint));
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp)
    grid on;
    zoom on;
    title('Cy5 Cy7 FRET');
    ylabel('FRET');
    xlabel('Time (s)');
    set(gca, 'FontSize', 12);
    
    subplot('position',[0.93 0.04 0.03 0.1]);
    x = -0.1:0.05:1.1;
    [hX,hN]=hist(trc.fret12(firstpoint:lastpoint),x);
    barh(hN,hX,'k');
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    xlabel('Count');
    grid on;
    axis on;
    zoom on;
end