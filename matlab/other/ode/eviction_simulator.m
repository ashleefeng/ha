figure
hold on
for tau1 = 60:60:120
    %tau1 = 120;
    for tau2 = tau1:tau1
    %tau2 = tau1;
    n = 10000;
    label_eff = 0.7;
    dead_frac = 0.2;
    time_first_evicted = zeros(n, 1);
    time_second_evicted = zeros(n, 1);
    num_cy3 = zeros(n, 2);
    dead_nuc = zeros(n, 1);
    
    count = 0;
    for i = 1:n
        
        rand1 = rand;
        rand2 = rand;
        rand3 = rand;
        if rand1 < label_eff
            num_cy3(i, 1) = 1; % label proximal H2A
        end
        if rand2 < label_eff
            num_cy3(i, 2) = 1; % label distal H2A
        end
        if num_cy3(i, 1) == 0 && num_cy3(i, 2) == 0
            count = count + 1;
        end
        if rand3 < dead_frac
            dead_nuc(i) = 1;
        end
    end
    percent_dark_init = 1 - count/n;
    
    for i = 1:n
        time_first_evicted(i) = exprnd(tau1);
        time_second_evicted(i) = exprnd(tau2) + time_first_evicted(i);
    end
    
    total_time = 3600;
    percent_fret = zeros(total_time + 1, 1);
    
    for t = 0:total_time
        dark_count = 0;
        for i = 1:n
            if num_cy3(i, 1) == 0 && num_cy3(i, 2) == 0
                dark_count = dark_count + 1;
            elseif time_second_evicted(i) < t && dead_nuc(i) ~= 1 % second eviction event
                dark_count = dark_count + 1;
            elseif time_first_evicted(i) < t && dead_nuc(i) ~= 1% first eviction event
                if rand < 0.5 % evicts proximal dimer
                    if num_cy3(i, 1) == 1 && num_cy3(i, 2) == 0
                        dark_count = dark_count + 1;
                    end
                else % evicts distal dimer
                    if num_cy3(i, 1) == 0 && num_cy3(i, 2) == 1
                        dark_count = dark_count + 1;
                    end
                end
            end
        end
        percent_fret(t+1) = (1 - dark_count / n) / percent_dark_init;
    end
    
    plot(0:1/60:60, percent_fret);
    end
end
xlabel('Time(min)');
ylabel('Fraction of nucleosomes with Cy3-H2A');
ylim([0, 1]);
legend('\tau = 60 s', '\tau = 120 s'); 

real = [[0, 1];[3, 0.54];[7, 0.43]; [11, 0.41]; [20, 0.34]; [30, 0.3]; [45, 0.29]; [60, 0.23]];
real2 = [[0, 1];[3, 0.63];[7, 0.49];[11, 0.4];[20, 0.33];[30, 0.38];[45, 0.27];[60, 0.24]];
plot(real(:, 1), real(:, 2), '.', 'MarkerSize', 10);
plot(real2(:,1), real2(:, 2), '.', 'MarkerSize', 10);
set(gca, 'FontSize', 12);