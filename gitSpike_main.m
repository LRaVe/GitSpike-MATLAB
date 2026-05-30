close all

% Make sure subfolders containing functions and order helpers are available
thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    repoRoot = fileparts(thisFile);
    addpath(genpath(repoRoot));
end


% ==== Selection showing + plotting ====
measures=63;               % +1:ISI,+2:SPIKE,+4:RI-SPIKE,+8:SPIKE-Synchro,+16:SPIKE-order,+32:Spike Train Order
adaptive_measures=15;       % +1:ISI,+2:SPIKE,+4:RI-SPIKE,+8:SPIKE-Synchro     % Adaptive
showing=14;                 % +1:Spike Trains,+2:Distance,+4:Profile,+8:Matrix
plotting=14;               % +1:Spike Trains,+2:Distance,+4:Profile,+8:Matrix
sort_spike_trains=0;       % 0-no,1-yes


% ==== Parameters ====
tmin=0;
tmax=10;
threshold=1000;

% ==== Dataset ====
dataset=4;

if dataset==3
    tmax=1;
    num_trains=3;
    spikes=cell(1,num_trains);
    spikes{1} = [0.0001 0.7142];
    spikes{2} = [0.2858 0.9999];                   
    spikes{3} = [0.1429 0.8571];
elseif dataset==4
    tmax=10;
    num_trains=3;
    spikes=cell(1,num_trains);
    spikes{1}=[0 1.9 3.9 7 10];
    spikes{2}=[0 2 7.1 9 10];
    spikes{3}=[0 2.1 4.1 6.9 10];
    %spikes{4}=[0 2.2 6.8 7.1 10];
end

plot_spikes = spikes;

% ==== Add auxiliary spikes at boundaries ====
number_spikes=sum(cellfun(@length,spikes));
[spikes, aux_begin, aux_end] = add_auxiliary_spikes(spikes, tmin, tmax);


% ==== SPIKE trains ====
if mod(showing,2)>0 
    for i=1:num_trains
        fprintf('Spike train %d: ', i);
        disp(spikes{i});
    end
end
if mod(plotting,2)>0
    figure;
    hold on;
    set(gcf, 'Color','w');
    for i=1:num_trains
        y = num_trains - i + 1;
        for j=1:numel(plot_spikes{i})
            line([plot_spikes{i}(j) plot_spikes{i}(j)], [y-0.35 y+0.35], 'Color', 'k', 'LineWidth', 1.5);
        end
    end
    set(gca, 'Color', 'w');
    set(gca, 'XColor', 'k');
    set(gca, 'YColor', 'k');
    set(gca, 'YTick', []);
    xMargin = 0.001;
    xlim([tmin - xMargin, tmax + xMargin]);
    ylim([0.5 num_trains+0.5]);
    hold off;
end

% ==== ISI distance ====
if mod(measures,2)>0
    f_ISI_distance(spikes, tmin, tmax, showing, plotting);
end

if mod(adaptive_measures,2)>0
    f_ISI_distance_adaptive_v1(spikes, tmin, tmax, threshold, showing, plotting);
end


% ==== SPIKE distance + RI-SPIKE ====
if mod(measures,8)>1 || mod(adaptive_measures,8)>1
    SPIKE_distances_all(spikes, tmin, tmax, threshold, measures, adaptive_measures, showing, plotting, aux_begin, aux_end);
end


% ==== SPIKE-Synchro ====
if mod(measures,16)>7 || mod(adaptive_measures,16)>7
    if mod(measures,16)>7
        [C_matrix, C_global, spike_synchro_times, spike_synchro_values] = f_spike_synchro_multi(spikes, tmin, tmax);
        spike_synchro_profile = [spike_synchro_times(:), spike_synchro_values(:)];
        if mod(showing,4)>1 || mod(plotting,16)>1
            if mod(showing,4)>1
                fprintf('\n=== Global SPIKE-Synchronization Index ===\n');
                fprintf('C_global: %.4f\n', C_global);
            end
        end
        if mod(showing,8)>3 || mod(plotting,8)>3
            if mod(showing,8)>3
                disp('=== SPIKE-Synchronization Profile ===');
                disp(spike_synchro_profile);
            end
                if mod(plotting,8)>3 && ~isempty(spike_synchro_profile)
                titleStr = sprintf('SPIKE-Synchronization C = %.4g', C_global);
                figure('Name', titleStr, 'NumberTitle', 'off');
                hold on;
                grid on;
                box on;
                plot([tmin,tmax],[C_global,C_global], '-', 'Color', 'red', 'LineWidth', 1);
                plot(spike_synchro_profile(:,1), spike_synchro_profile(:,2), '-o', 'Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 6);
                xlim([tmin,tmax]);
                ylim([0,1]);
                title(titleStr);
                hold off;
            end
        end
        if mod(showing,16)>7 || mod(plotting,16)>7
            if mod(showing,16)>7
                disp('=== Pairwise Coincidence Matrix ===');
                disp(C_matrix);
            end
                if mod(plotting,16)>7
                titleStr = sprintf('Pairwise Coincidence Matrix (C_{global} = %.4f)', C_global);
                figure('Name', titleStr, 'NumberTitle', 'off');
                n = length(spikes);
                matrix_min = min(C_matrix(:));
                matrix_max = max(C_matrix(:));
                if matrix_min == matrix_max
                    matrix_max = matrix_min + 1;
                end
                imagesc(C_matrix, [matrix_min matrix_max]);
                colormap(jet);
                colorbar;
                box on;
                axis equal;
                xlim([0.5 n+0.5]);
                ylim([0.5 n+0.5]);
                set(gca, 'XDir', 'normal');
                set(gca, 'YDir', 'reverse');
                set(gca, 'XTick', 1:n, 'YTick', 1:n);
                xlabel('Spike trains');
                ylabel('Spike trains');
                title(titleStr);
            end
        end
    end

    if mod(adaptive_measures,16)>7
        [C_matrix_adaptive, C_global_adaptive, spike_synchro_times_adaptive, spike_synchro_values_adaptive] = f_spike_synchro_multi(spikes, tmin, tmax, 'auto');
        spike_synchro_profile_adaptive = [spike_synchro_times_adaptive(:), spike_synchro_values_adaptive(:)];
        if mod(showing,4)>1 || mod(plotting,16)>1
            if mod(showing,4)>1
                fprintf('\n=== Global Adaptive SPIKE-Synchronization Index ===\n');
                fprintf('C_global: %.4f\n', C_global_adaptive);
            end
        end
        if mod(showing,8)>3 || mod(plotting,8)>3
            if mod(showing,8)>3
                disp('=== Adaptive SPIKE-Synchronization Profile ===');
                disp(spike_synchro_profile_adaptive);
            end
                if mod(plotting,8)>3 && ~isempty(spike_synchro_profile_adaptive)
                titleStr = sprintf('Adaptive SPIKE-Synchronization C = %.4g', C_global_adaptive);
                figure('Name', titleStr, 'NumberTitle', 'off');
                hold on;
                grid on;
                box on;
                plot([tmin,tmax],[C_global_adaptive,C_global_adaptive], '-', 'Color', 'red', 'LineWidth', 1);
                plot(spike_synchro_profile_adaptive(:,1), spike_synchro_profile_adaptive(:,2), '-o', 'Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 6);
                xlim([tmin,tmax]);
                ylim([0,1]);
                title(titleStr);
                hold off;
            end
        end
        if mod(showing,16)>7 || mod(plotting,16)>7
            if mod(showing,16)>7
                disp('=== Pairwise Adaptive Coincidence Matrix ===');
                disp(C_matrix_adaptive);
            end
                if mod(plotting,16)>7
                titleStr = sprintf('Pairwise Adaptive Coincidence Matrix (C_{global} = %.4f)', C_global_adaptive);
                figure('Name', titleStr, 'NumberTitle', 'off');
                n = length(spikes);
                matrix_min = min(C_matrix_adaptive(:));
                matrix_max = max(C_matrix_adaptive(:));
                if matrix_min == matrix_max
                    matrix_max = matrix_min + 1;
                end
                imagesc(C_matrix_adaptive, [matrix_min matrix_max]);
                colormap(jet);
                colorbar;
                box on;
                axis equal;
                xlim([0.5 n+0.5]);
                ylim([0.5 n+0.5]);
                set(gca, 'XDir', 'normal');
                set(gca, 'YDir', 'reverse');
                set(gca, 'XTick', 1:n, 'YTick', 1:n);
                xlabel('Spike trains');
                ylabel('Spike trains');
                title(titleStr);
            end
        end
    end
end


% ==== SPIKE-Order ====
if mod(measures,32)>15
    [sortedOrders,sortedTimes,SO_matrix]=order_spikes(tmin,tmax,spikes);
    if mod(showing,4)>1 || mod(plotting,16)>1
        spike_order=0;
        if mod(showing,4)>1
            disp(['Spike order D = ', num2str(spike_order)]);
        end
    end
    if mod(showing,8)>3 || mod(plotting,8)>3 
        if mod(showing,8)>3 % showing profile
            spike_order_profile = [sortedTimes(:), sortedOrders(:)];
            disp('Spike order profile: ');
            disp(spike_order_profile);
        end
        if mod(plotting,8)>3 % plotting profile
            titleStr = sprintf('Spike order D = %.4g', spike_order);
            figure('Name', titleStr, 'NumberTitle', 'off');
            hold on;
            grid on;
            box on;
            plot([tmin,tmax],[spike_order,spike_order], '-', 'Color', 'red', 'LineWidth', 1);
            plot(sortedTimes,sortedOrders,'-o','Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 6);
            xlim([tmin,tmax]);
            ylim([-1.1,1.1]);
            title(titleStr);
            yticks([-1,0,1]);
            hold off;
        end
    end
    if mod(showing,16)>7 || mod(plotting,16)>7
        if mod(showing,16)>7 % showing matrix
            disp('Spike order matrix: ');
            disp(SO_matrix);
        end
        if mod(plotting,16)>7 % plotting matrix
            n=length(spikes);
            titleStr = sprintf('Spike order matrix D = %g', spike_order);
            figure('Name', titleStr, 'NumberTitle', 'off');
            hold on;
            imagesc(SO_matrix, [-1 1]);
            box on;
            colormap(jet);
            colorbar;
            axis equal;
            xlim([0.5 n+0.5]);
            ylim([0.5 n+0.5]);
            set(gca, 'XDir', 'normal');
            set(gca, 'YDir', 'reverse');
            set(gca, 'XTick', 1:n, 'YTick', 1:n);
            xlabel('Spike trains');
            ylabel('Spike trains');
            title(titleStr);
            hold off;
        end
    end
end

% ==== SPIKE-Train-Order ====
if mod(measures,64)>31
    [results,order_matrix]=order_trains(tmin,tmax,spikes);
    if mod(showing,16)>1 || mod(plotting,16)>1
        [F,sortedTimes,sortedOrders]=compute_spike_train_order_value(spikes,results,number_spikes);
        if mod(showing,4)>1
            fprintf('Spike train order F = %.4f\n', F);
        end
    end
    if mod(showing,8)>3 || mod(plotting,8)>3
        if mod(showing,8)>3 % showing profile
            spike_train_order_profile = [sortedTimes(:), sortedOrders(:)];
            disp('Spike train order profile: ');
            disp(spike_train_order_profile);
        end
        if mod(plotting,8)>3 % plotting profile
            titleStr = sprintf('Spike train order F = %.4g', F);
            figure('Name', titleStr, 'NumberTitle', 'off');
            hold on;
            grid on;
            box on;
            plot([tmin,tmax],[F,F], '-', 'Color', 'red', 'LineWidth', 1);
            plot(sortedTimes,sortedOrders,'-o','Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 6);
            xlim([tmin,tmax]);
            ylim([-1.1,1.1]);
            title(titleStr);
            yticks([-1,0,1]);
            hold off;
        end
    end
    if mod(showing,16)>7 || mod(plotting,16)>7
        if mod(showing,16)>7 % showing matrix
            disp('Spike train order matrix: ');
            disp(order_matrix);
        end
        if mod(plotting,16)>7 % plotting matrix
            n=length(spikes);
            titleStr = sprintf('Spike train order F = %g', F);
            figure('Name', titleStr, 'NumberTitle', 'off');
            hold on;
            matrix_min = min(order_matrix(:));
            matrix_max = max(order_matrix(:));
            imagesc(order_matrix, [matrix_min matrix_max]);
            box on;
            colormap(jet);
            colorbar;
            axis equal;
            xlim([0.5 n+0.5]); 
            ylim([0.5 n+0.5]);
            set(gca, 'XDir', 'normal');
            set(gca, 'YDir', 'reverse');
            set(gca, 'XTick', 1:n, 'YTick', 1:n);
            xlabel('Spike trains');
            ylabel('Spike trains');
            title(titleStr);
            hold off;
        end
    end
end 


    
