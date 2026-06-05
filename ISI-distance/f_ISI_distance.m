%% ISI-distance computation with auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026
function f_ISI_distance(spikes, tmin, ...
    tmax, showing, plotting)
% F_ISI_DISTANCE Calculates the ISI-distance between multiple spike trains
%
% Inputs:
%   spikes_trains : cell array containing spike timing vectors
%   tmin, tmax    : time boundaries (minimum and maximum)
%   showing       : Control console output (+2: Pairs, +4: Pop Profile, +8: Matrix)
%   plotting      : Control figures (+2: Pairs, +4: Pop Profile, +8: Matrix)

    % Manage parameters if necessary 
    if nargin < 4 || isempty(showing), showing = 15; end
    if nargin < 5 || isempty(plotting), plotting = 15; end
    num_trains = length(spikes);
    
    if num_trains < 2
        if bitand(showing, 2)
            disp('Not enough spike trains to calculate a distance.');
        end
        return;
    end
    
    num_pairs = (num_trains * (num_trains - 1)) / 2;
    dist_matrix = zeros(num_trains, num_trains);
    I = zeros(1, num_pairs); 
    pair_data = cell(1, num_pairs); 
    
    % % Edge correction 
    % spikes = cell(1, num_trains);
    % for i = 1:num_trains
    %     s = unique(spikes_trains{i}); 
    %     spikes{i} = s(s >= tmin & s <= tmax); 
    % end
    
    compteur = 0;
    num_cols = 2; 
    num_rows = ceil(num_pairs / num_cols);
    
    % Set figure window name for pairwise plots
    if bitand(plotting, 4)
        figure('Name', 'Pairwise ISI Distances (Classic Mode)');
        set(gcf, 'Name', 'Pairwise ISI Distances (Classic Mode)');
    end
    
    for i = 1:num_trains
        for j = i+1:num_trains
            compteur = compteur + 1;
            t_all = unique([spikes{i}, spikes{j}]);
            Iij = 0;
            It_list = zeros(1, length(t_all)-1); 

            for k = 1 : length(t_all)-1
                t_mid = (t_all(k) + t_all(k+1)) / 2;
                
                % Train i
                idx = find(spikes{i} <= t_mid, 1, 'last');
                if isempty(idx), idx = 1; end
                if idx >= length(spikes{i}), idx = length(spikes{i}) - 1; end
                
                val_x = spikes{i}(idx+1) - spikes{i}(idx); 
                
                % Train j
                idy = find(spikes{j} <= t_mid, 1, 'last');
                if isempty(idy), idy = 1; end
                if idy >= length(spikes{j}), idy = length(spikes{j}) - 1; end
                
                val_y = spikes{j}(idy+1) - spikes{j}(idy);
               
                if isempty(val_x) || val_x < 0, val_x = 0; end
                if isempty(val_y) || val_y < 0, val_y = 0; end
                
                % Calcul of the ISI distance (avoid the division by 0)
                if max(val_x, val_y) > 0
                    I_t = abs(val_x - val_y) / max(val_x, val_y);
                else
                    I_t = 0;
                end
                It_list(k) = I_t;
                
                % Integration between the realtime window
                segment_tmin = max(t_all(k), tmin);
                segment_tmax = min(t_all(k+1), tmax);
                if segment_tmax > segment_tmin
                    Iij = Iij + I_t * (segment_tmax - segment_tmin);
                end
            end

            % for k = 1 : length(t_all)-1
            %     t_mid = (t_all(k) + t_all(k+1)) / 2;
            % 
            %     % %% Edge correction %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     % % train i
            %     % if isempty(spikes{i})
            %     %     val_x = tmax - tmin; 
            %     % elseif t_mid < spikes{i}(1)
            %     %     val_x = spikes{i}(1) - tmin; 
            %     % elseif t_mid > spikes{i}(end)
            %     %     val_x = tmax - spikes{i}(end); 
            %     % else
            %     %     idx = find(spikes{i} <= t_mid, 1, 'last');
            %     %     val_x = spikes{i}(idx+1) - spikes{i}(idx); 
            %     % end
            %     % 
            %     % % train j
            %     % if isempty(spikes{j})
            %     %     val_y = tmax - tmin; 
            %     % elseif t_mid < spikes{j}(1)
            %     %     val_y = spikes{j}(1) - tmin;
            %     % elseif t_mid > spikes{j}(end)
            %     %     val_y = tmax - spikes{j}(end);
            %     % else
            %     %     idy = find(spikes{j} <= t_mid, 1, 'last');
            %     %     val_y = spikes{j}(idy+1) - spikes{j}(idy);
            %     % end
            %     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 
            %     % Train i
            %     idx = find(spikes{i} <= t_mid, 1, 'last');
            %     val_x = spikes{i}(idx+1) - spikes{i}(idx); 
            % 
            %     % Train j
            %     idy = find(spikes{j} <= t_mid, 1, 'last');
            %     val_y = spikes{j}(idy+1) - spikes{j}(idy);
            % 
            %     I_t = abs(val_x - val_y) / max(val_x, val_y);
            %     Iij = Iij + I_t * (t_all(k+1) - t_all(k));
            %     It_list(k) = I_t;
            % end
            
            I_final_pair = Iij / (tmax - tmin);
            dist_matrix(i,j) = I_final_pair;
            dist_matrix(j,i) = I_final_pair;
            I(compteur) = I_final_pair; 
            
            pair_data{compteur}.t = t_all;
            pair_data{compteur}.It = It_list;
            
            if bitand(plotting, 4) 
                subplot(num_rows, num_cols, compteur); 
                I_plot = [It_list, It_list(end)];
                stairs(t_all, I_plot, 'LineWidth', 1.5); 
                title(['Pair ', num2str(i), ' & ', num2str(j)]);  
                subtitle(['Dist: ', num2str(I_final_pair, '%.4f')]); 
                xlabel('Time'); 
                ylabel('I(t)'); 
                xlim([tmin tmax]); 
                ylim([0 1]); 
                box on; grid on; 
            end
        end
    end   
    
    I_mean = mean(I);
    if bitand(showing, 2)
        fprintf('The ISI-distance is: %.4f\n', I_mean);
    end
 
    all_spikes_combined = [spikes{:}]; 
    %t_global = unique([tmin, all_spikes_combined, tmax]);
    t_global = unique(all_spikes_combined);
    
    % Population profile average matrix allocation
    I_matrix = zeros(length(pair_data), length(t_global)-1);
    for p = 1:length(pair_data)
        t_p = pair_data{p}.t;
        It_p = pair_data{p}.It;
        for k = 1:length(t_global)-1
            t_mid = (t_global(k) + t_global(k+1)) / 2;
            idx = find(t_p(1:end-1) <= t_mid, 1, 'last');
            if isempty(idx)
                idx = 1;
            end
            I_matrix(p, k) = It_p(idx);
        end
    end
    I_pop_mean = mean(I_matrix, 1);
    
    % Matrix plot 
    if bitand(plotting, 8)
        title_mat = sprintf('Matrix of the ISI-distance - Population Mean: %.4f', I_mean);
        figure('Name', title_mat);
        set(gcf, 'Name', title_mat);
        imagesc(dist_matrix); 
        colorbar;
        colormap jet;
        title(title_mat);
        xlabel('Spike Train Index'); ylabel('Spike Train Index');
        box on;
    end
    
    if bitand(showing, 8)
        disp('Final ISI-Distance matrix:');
        disp(dist_matrix);
    end
    
    % Population Profile plot 
    if bitand(plotting, 4)
        title_pop = sprintf('Evolution of Population Average ISI distance - Global: %.4f', I_mean);
        figure('Name', title_pop);
        set(gcf, 'Name', title_pop); 
        stairs(t_global, [I_pop_mean, I_pop_mean(end)], 'LineWidth', 1.5);
        xlabel('Time'); 
        ylabel('Average I(t)');
        xlim([0 tmax]);   
        ylim([0 1]);
        title(title_pop);
        box on; 
        grid on;
    end
    
    if bitand(showing, 4)
        fprintf('\n=== Final ISI-Distance plot : ===\n');
        fprintf('  Time(t)  |  Average ISI Distance I(t)\n');
        I_pop_extended = [I_pop_mean, I_pop_mean(end)];
        for idx_plot = 1:length(t_global)
            fprintf('      %8.4f     |      %8.4f\n', t_global(idx_plot), ...
                I_pop_extended(idx_plot));
        end
    end
end