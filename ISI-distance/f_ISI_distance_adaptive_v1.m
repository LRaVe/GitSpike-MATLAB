%% ISI-distance adaptive computationwith auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026


function f_ISI_distance_adaptive_v1(spikes, ...
    tmin, tmax, threshold, showing, plotting)
% COMPUTE_ADAPTIVE_ISI Calculates the Adaptive ISI-distance between spike trains
%
% Inputs:
%   spikes_trains : cell array containing spike timing vectors
%   tmin, tmax    : time boundaries (minimum and maximum)
%   threshold     : 0 for classic, >0 for manual MRTS, 'auto' for automated MRTS
%   showing       : Control console output (+2: Pairs, +4: Pop Profile, +8: Matrix)
%   plotting      : Control figures (+2: Pairs, +4: Pop Profile, +8: Matrix)

    % Manage of the parameters if necessary
    if nargin < 4 || isempty(threshold), threshold = 0; end
    if nargin < 5 || isempty(showing), showing = 15; end   
    if nargin < 6 || isempty(plotting), plotting = 15; end
    
    num_trains = length(spikes);
    if num_trains < 2
        if bitand(showing, 2)
            disp('Not enought train to calculate the ISI-distance');
        end
        return;
    end
    
    num_pairs = (num_trains * (num_trains - 1)) / 2;
    
    %% Dynamic figure number allocation
    fig_profile_id = 151;
    fig_matrix_id = 152;
    if bitand(plotting, 12)
        all_figs = findall(0, 'Type', 'figure'); 
        
        if ~isempty(all_figs)
           % Guarantees 'existing_ids' exists and is initialized as empty if 
           % no valid figures exist
            existing_ids = [all_figs.Number];
        else
            existing_ids = [];
        end
        while any(existing_ids == fig_profile_id) || any(existing_ids == fig_matrix_id)
                fig_profile_id = fig_profile_id + 2; 
                fig_matrix_id = fig_matrix_id + 2;  
        end
    end 

    % % Edge correction 
    % spikes = cell(1,num_trains);
    % for i = 1:num_trains
    %     s = unique(spikes_trains{i}); 
    %     spikes{i} = s(s >= tmin & s <= tmax); 
    % end
    
    if ischar(threshold) && strcmpi(threshold, 'auto')
        MRTS = autoMRTS(spikes);
    else
        MRTS = threshold; 
    end
    
    if MRTS > 0
        mode_label = sprintf('Adaptive (MRTS = %.3f)', MRTS);
    else
        mode_label = 'Classic (MRTS = 0)';
    end
    
    dist_matrix = zeros(num_trains, num_trains);
    I = zeros(1, num_pairs);
    pair_data = cell(1, num_pairs); 
    
    compteur = 0;
    %num_cols = 2; 
    %num_rows = ceil(num_pairs / num_cols);
    
    % if bitand(plotting, 4)
    %     title_pairs = ['Pairwise ISI Distances - ' mode_label];
    %     figure('Name', title_pairs);
    %     set(gcf, 'Name', title_pairs);
    % end
    
    for i = 1:num_trains
        for j = i+1:num_trains
            compteur = compteur + 1;
            t_all = unique([spikes{i}, spikes{j}]);
            It_list = zeros(1, length(t_all)-1);

            for k = 1 : length(t_all)-1
                t_mid = (t_all(k) + t_all(k+1)) / 2;
                
                % Train i 
                idx = find(spikes{i} <= t_mid, 1, 'last');
                if isempty(idx), idx = 1; end
                if idx >= length(spikes{i}), idx = length(spikes{i}) - 1; end
                
                vx = spikes{i}(idx+1) - spikes{i}(idx); 
                
                % Train j 
                idy = find(spikes{j} <= t_mid, 1, 'last');
                if isempty(idy), idy = 1; end
                if idy >= length(spikes{j}), idy = length(spikes{j}) - 1; end
                
                vy = spikes{j}(idy+1) - spikes{j}(idy);
                
                % Avoid nul values and division by 0
                if isempty(vx) || vx < 0, vx = 0; end
                if isempty(vy) || vy < 0, vy = 0; end
                
                denominateur = max([vx, vy, MRTS]);
                if denominateur > 0
                    It_list(k) = abs(vx - vy) / denominateur;
                else
                    It_list(k) = 0;
                end
            end

            Iij = 0;
            for k = 1:length(t_all)-1
                segment_tmin = max(t_all(k), tmin);
                segment_tmax = min(t_all(k+1), tmax);
                
                if segment_tmax > segment_tmin
                    Iij = Iij + It_list(k) * (segment_tmax - segment_tmin);
                end
            end
            
            Iij = Iij / (tmax - tmin);
            dist_matrix(i,j) = Iij;
            dist_matrix(j,i) = Iij;
            I(compteur) = Iij;
            
            pair_data{compteur}.t = t_all;
            pair_data{compteur}.It = It_list;
        end
    end   
    
    I_mean = mean(I);
    if bitand(showing, 2)
        fprintf('The global ISI-distance is: %.4f\n', I_mean);
    end
    
    % Average of the population
    all_spikes_combined = [spikes{:}]; 
    t_global = unique(all_spikes_combined); 
    I_matrix = zeros(length(pair_data), length(t_global)-1);
    for p = 1:length(pair_data)
        for k = 1:length(t_global)-1
            t_mid = (t_global(k) + t_global(k+1)) / 2;
            idx = find(pair_data{p}.t(1:end-1) <= t_mid, 1, 'last');
            if isempty(idx)
                idx = 1;
            end
            I_matrix(p, k) = pair_data{p}.It(idx);
        end
    end
    I_pop_mean = mean(I_matrix, 1);

    if bitand(plotting, 4)
        title_pop = ['Population Average - ' mode_label ' Population Mean: ' num2str(I_mean, '%.4f')];
        figure(fig_profile_id);
        set(gcf, 'Name', title_pop); 
        stairs(t_global, [I_pop_mean, I_pop_mean(end)], 'LineWidth', 1.5);
        title('Population Average');
        subtitle(['Global ISI-distance: ', num2str(I_mean, '%.4f')]);
        xlabel('Time');
        ylabel('Average I(t)'); 
        xlim([0 tmax]); 
        ylim([0 1]); 
        box on; 
        grid on;
    end
       
    if bitand(showing, 4)
        fprintf('\n=== Final ISI-Distance plot (%s) : ===\n', mode_label);
        fprintf('  Time(t)  |  Average ISI Distance I(t)\n');
        I_pop_extended = [I_pop_mean, I_pop_mean(end)];
        for idx_plot = 1:length(t_global)
            fprintf('      %8.4f     |      %8.4f\n', t_global(idx_plot), ...
                I_pop_extended(idx_plot));
        end
    end
    
    % Global plots
    if bitand(plotting, 8)
        title_mat = ['ISI Matrix - ' mode_label ' Population Mean: ' num2str(I_mean, '%.4f')];
        figure(fig_matrix_id);
        set(gcf, 'Name', title_mat); 
        imagesc(dist_matrix); 
        colorbar;
        colormap jet; 
        title('ISI Matrix');
        subtitle(['Global ISI-distance: ', num2str(I_mean, '%.4f')]);
        xlabel('Spike Train Index'); 
        ylabel('Spike Train Index');
        box on; 
    end
    
    if bitand(showing, 8)
        disp('Final ISI-Distance matrix:');
        disp(dist_matrix);
    end
end

% Subfunction for the MRTS parameter
function [MRTS] = autoMRTS(spikes)
    sum_isi_sqr = 0;
    num_isi = 0;
    for i=1:length(spikes)
        if length(spikes{i}) >= 2
            for j=1:(length(spikes{i})-1)
                sum_isi_sqr = sum_isi_sqr + (spikes{i}(j+1)-spikes{i}(j))^2;
                num_isi = num_isi + 1;
            end
        end
    end
    if num_isi > 0
        MRTS = (sum_isi_sqr/num_isi)^0.5;
    else
        MRTS = 0;
    end
end