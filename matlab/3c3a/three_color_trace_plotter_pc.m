function three_color_trace_plotter_pc(trc, trc_name, firstpoint, lastpoint)
    
    ex0_ylim = [-100, 500];
    ex1_ylim = [-100, 700];
    ex2_ylim = [-50, 300];
    fret01_end = floor(500*300/300);
    fret02_end = floor(500*300/300);
    %fret01_end = lastpoint;
    %fret02_end = lastpoint;
    fret12_end = lastpoint;
    
    %% green laser excitation corrected trace
    
    subplot('position',[0.1 0.84 0.8 0.14]);
    set(gcf,'Position',[400 400 400 500]); % standard traces window
    %set(gcf, 'Position', [400 400, 1000 500]); % expanded
    %set(gcf, 'Position', [400 400 200 500]); % shrinked
    
    plot(trc.time(firstpoint:lastpoint), trc.ex0em0(firstpoint:lastpoint), 'g', trc.time(firstpoint:lastpoint), trc.ex0em1(firstpoint:lastpoint), 'r', 'LineWidth', 1);
    hold on
    plot(trc.time(firstpoint:lastpoint), trc.ex0em2(firstpoint:lastpoint), 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1);
    %plot(trc.time, trc.ex0em0 + trc.ex0em1 + trc.ex0em2 + 200, 'k');
    hold off
    %title(['Green Laser Molecule ' trc_name], 'Interpreter', 'none');
    %title('Cy3 Excitation');
    %ylabel('Intensity');
    zoom on;
    
    set(gca, 'FontSize', 8);    
    ylim(ex0_ylim);
    %xlim([130 230]);
%     xticks([0 20 40 60 80 100]);
    set(gca,'XTickLabel',[]);
    
    %% red laser excitation corrected trace
    
    subplot('position',[0.1 0.68 0.8 0.14]); 
    %plot(trc.time, trc.ex1em0, 'g', trc.time, trc.ex1em1, 'r', trc.time, trc.ex1em2, 'm');
    plot(trc.time(firstpoint:lastpoint), trc.ex1em1(firstpoint:lastpoint), 'r', 'LineWidth', 1);
    hold on
    plot(trc.time(firstpoint:lastpoint), trc.ex1em2(firstpoint:lastpoint), 'Color',[0.4940, 0.1840, 0.5560], 'LineWidth', 1);
    hold off
    %title(['Red Laser Molecule ' trc_name], 'Interpreter', 'none');
    %title('Cy5 Excitation');
    zoom on;
    %ylabel('Intensity');
    
    set(gca, 'FontSize', 8);
    ylim(ex1_ylim);
   % xlim([130 230]);
%     xticks([0 20 40 60 80 100]);
    set(gca,'XTickLabel',[]);
    %yticks([0 400 800]);
    %set(gca,'YTickLabel',['', 0, 200, 400, 600]);
    
    %% 750 laser excitation corrected trace
    
    subplot('position',[0.1 0.52 0.8 0.14]);
    %plot(trc.time, trc.ex2em0, 'g', trc.time, trc.ex2em1, 'r', trc.time, trc.ex2em2, 'm');
    plot(trc.time(firstpoint:lastpoint), trc.ex2em2(firstpoint:lastpoint), 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1);
    %title(['750 Laser Molecule ' trc_name], 'Interpreter', 'none');
    %title('Cy7 Excitation');
    zoom on;
    %ylabel('Intensity');
    
    set(gca, 'FontSize', 8);
    ylim(ex2_ylim);
    %xlim([130 230]);
%     xticks([0 20 40 60 80 100]);
    set(gca,'XTickLabel',[]);

%% cy3-cy5 fret trace

    subplot('position',[0.1 0.36 0.8 0.14]);
    plot(trc.time(firstpoint:fret01_end), trc.fret01(firstpoint:fret01_end), 'LineWidth', 1);
    ylim([-0.1 1.1]);
    
%     xlim([0 100]);
%     xticks([0 20 40 60 80 100]);
    set(gca, 'FontSize', 8);
    set(gca,'XTickLabel',[]);
    
%% cy3-cy7 fret trace

    subplot('position',[0.1 0.20 0.8 0.14]);
    plot(trc.time(firstpoint:fret02_end), trc.fret02(firstpoint:fret02_end), 'LineWidth', 1);
    ylim([-0.1 1.1]);
    %grid on;
    zoom on;
    
%     xlim([0 100]);
%     xticks([0 20 40 60 80 100]);
    set(gca, 'FontSize', 8);
    set(gca,'XTickLabel',[]);

%% cy5-cy7 fret trace

    subplot('position', [0.1 0.04 0.8 0.14]);
    plot(trc.time(firstpoint:fret12_end), trc.fret12(firstpoint:fret12_end), 'LineWidth', 1);
    %hold on
    %plot(trc.time(firstpoint:fret02_end), trc.fret02(firstpoint:fret02_end), 'LineWidth', 1);
    %hold off
    ylim([-0.1 1.1]);
    %grid on;
    zoom on;
    
    %title('Cy5 Cy7 FRET');
    %ylabel('FRET');
    xlabel('Time (s)');
    set(gca, 'FontSize', 8);
    %set(gca,'XTickLabel',[]);
    %xlim([130 230]);
%     xticks([0 20 40 60 80 100]);
end