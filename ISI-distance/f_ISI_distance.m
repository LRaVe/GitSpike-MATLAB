%% ISI-distance computation with auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026

function f_ISI_distance(spikes, tmin,tmax, showing, plotting)
% F_ISI_DISTANCE Computes the bivariate and multivariate ISI-distance between spike trains.
%
%   Calculates the instantaneous and overall Inter-Spike Interval (ISI) distance 
%   between multiple spike trains within a given temporal window.
%
%   .. note::
%      Based on the ISI-distance algorithm proposed by Kreuz et al.
%
%   .. seealso::
%      Kreuz T, et al. *Measuring spike train synchrony.* J Neurosci Methods, 2007.
%
%   The instantaneous bivariate ISI-distance profile :math:`I(t)` between two spike trains is defined as:
%
%   .. math::
%
%      I(t) = \frac{|x_1(t) - x_2(t)|}{\max(x_1(t), x_2(t))}
%
%   where :math:`x_1(t)` and :math:`x_2(t)` are the instantaneous inter-spike intervals 
%   for each spike train at time :math:`t`.
%
%   The overall bivariate ISI-distance :math:`D_I` is the time average over the window :math:`[T_0, T_1]`:
%
%   .. math::
%
%      D_I = \frac{1}{T_1 - T_0} \int_{T_0}^{T_1} I(t) \, dt
%
%   In the multivariate case, the population profile :math:`I_{multi}(t)` is averaged over all :math:`N(N-1)/2` pairs:
%
%   .. math::
%
%      I_{multi}(t) = \frac{2}{N(N-1)} \sum_{<i,j>} I^{i,j}(t)
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Basic call with default output flags (showing=15, plotting=15)
%      f_ISI_distance(spikes, tmin, tmax);
%
%      % Suppress all plots and console outputs
%      f_ISI_distance(spikes, tmin, tmax, 0, 0);
%
%      % Display only the population profile (console and plot)
%      f_ISI_distance(spikes, tmin, tmax, 4, 4);
%
%   :param spikes: Cell array where each element is a vector containing spike timestamps for a train.
%   :type spikes: cell
%   :param tmin: Lower temporal boundary of the window of interest.
%   :type tmin: double
%   :param tmax: Upper temporal boundary of the window of interest.
%   :type tmax: double
%   :param showing: Bitmask controlling console output verbosity (default: 15).
%                   * +2: Print pairwise ISI distance summary
%                   * +4: Print population profile array
%                   * +8: Print final pairwise distance matrix
%   :type showing: integer, optional
%   :param plotting: Bitmask controlling figure generation (default: 15).
%                    * +2: (Reserved for pairwise plots)
%                    * +4: Plot population profile over time
%                    * +8: Plot final ISI distance matrix heatmap
%   :type plotting: integer, optional
%
%   .. note::
%      The function requires at least two spike trains (`length(spikes) >= 2`) 
%      to compute pairwise distances.

    % Manage parameters if necessary 
    if nargin < 4 || isempty(showing), showing = 15; end
    if nargin < 5 || isempty(plotting), plotting = 15; end
    num_trains = length(spikes);

    %% Dynamic figure number allocation
    fig_profile_id = 101;
    fig_matrix_id = 102;
    if bitand(plotting, 12)
        all_figs = findall(0, 'Type', 'figure'); 
        % Guarantees 'existing_ids' exists and is initialized as empty if 
        % no valid figures exist
        if ~isempty(all_figs)
            existing_ids = [all_figs.Number];
        else
            existing_ids = [];
        end
        while any(existing_ids == fig_profile_id) || any(existing_ids == fig_matrix_id)
                fig_profile_id = fig_profile_id + 2; 
                fig_matrix_id = fig_matrix_id + 2;  
        end
    end
    

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
    compteur = 0;
    
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

            I_final_pair = Iij / (tmax - tmin);
            dist_matrix(i,j) = I_final_pair;
            dist_matrix(j,i) = I_final_pair;
            I(compteur) = I_final_pair; 
            
            pair_data{compteur}.t = t_all;
            pair_data{compteur}.It = It_list;
        end
    end   
    
    I_mean = mean(I);
    if bitand(showing, 2)
        fprintf('The ISI-distance is: %.4f\n', I_mean);
    end
 
    all_spikes_combined = [spikes{:}]; 
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

    % Population Profile plot 
    if bitand(plotting, 4)
        title_pop = sprintf('Evolution of Population Average ISI distance - Global: %.4f', I_mean);
        figure(fig_profile_id);
        set(gcf, 'Name', title_pop);
        stairs(t_global, [I_pop_mean, I_pop_mean(end)], 'LineWidth', 1.5);
        xlabel('Time'); 
        ylabel('Average I(t)');
        xlim([0 tmax]);   
        ylim([0 1]);
        title('Population Average');
        subtitle(['Global ISI-distance: ', num2str(I_mean, '%.4f')]);
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
    
    % Matrix plot 
    if bitand(plotting, 8)
        title_mat = sprintf('Matrix of the ISI-distance - Population Mean: %.4f', I_mean);
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