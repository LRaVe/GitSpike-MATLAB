%% ISI-distance adaptive computation with auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026

function f_ISI_distance_adaptive_v1(spikes, tmin, tmax, threshold, showing, plotting)
% F_ISI_DISTANCE_ADAPTIVE_V1 Computes the adaptive ISI-distance between spike trains.
%
%   Calculates the adaptive Inter-Spike Interval (ISI) distance between 
%   multiple spike trains using a Minimum Resolvable Time Scale (MRTS) parameter. 
%   The MRTS threshold prevents division-by-zero artifacts and stabilizes the distance 
%   measure during high firing rate periods or silent intervals.
%
%   The instantaneous adaptive ISI-distance profile :math:`I_{adapt}(t)` is defined as:
%
%   .. math::
%
%      I_{adapt}(t) = \frac{|x_1(t) - x_2(t)|}{\max\left(x_1(t), x_2(t), \tau_{MRTS}\right)}
%
%   where :math:`x_1(t)` and :math:`x_2(t)` are the instantaneous ISIs of the two spike trains 
%   at time :math:`t`, and :math:`\tau_{MRTS}` is the user-defined or automated threshold parameter.
%
%   The overall adaptive distance :math:`D_I` is obtained by time-averaging over :math:`[T_0, T_1]`:
%
%   .. math::
%
%      D_I = \frac{1}{T_1 - T_0} \int_{T_0}^{T_1} I_{adapt}(t) \, dt
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Classic mode (MRTS = 0)
%      f_ISI_distance_adaptive_v1(spikes, tmin, tmax);
%
%      % Automatic threshold selection based on RMS of ISIs
%      f_ISI_distance_adaptive_v1(spikes, tmin, tmax, 'auto');
%
%      % Manual MRTS threshold set to 0.05 seconds
%      f_ISI_distance_adaptive_v1(spikes, tmin, tmax, 0.05);
%
%      % Custom threshold with console/plot display bitmasks
%      f_ISI_distance_adaptive_v1(spikes, tmin, tmax, 'auto', 15, 15);
%
%   :param spikes: Cell array where each element is a vector containing spike timestamps for a train.
%   :type spikes: cell
%   :param tmin: Lower temporal boundary of the analysis window.
%   :type tmin: double
%   :param tmax: Upper temporal boundary of the analysis window.
%   :type tmax: double
%   :param threshold: MRTS threshold parameter (default: 0).
%                     * ``0``: Classic ISI-distance mode.
%                     * ``>0``: Manual MRTS scalar value.
%                     * ``'auto'``: Automated calculation based on the Root Mean Square (RMS) of all ISIs.
%   :type threshold: double or char, optional
%   :param showing: Bitmask controlling console output verbosity (default: 15).
%                   * +2: Print global ISI distance summary
%                   * +4: Print population profile array
%                   * +8: Print final distance matrix
%   :type showing: integer, optional
%   :param plotting: Bitmask controlling figure visualization (default: 15).
%                    * +2: (Reserved for pairwise plots)
%                    * +4: Plot population profile over time
%                    * +8: Plot distance matrix heatmap
%   :type plotting: integer, optional
%
%   .. note::
%      If `threshold` is set to ``'auto'``, the function internally calculates 
%      the Root Mean Square (RMS) of all inter-spike intervals:
%
%      .. math::
%
%         \tau_{MRTS} = \sqrt{\frac{1}{N_{isi}} \sum_{k=1}^{N_{isi}} (\Delta t_k)^2}
%
%   :Author: Laure WOLFF
%   :Date: May 2026

    % Manage the parameters if necessary
    if nargin < 4 || isempty(threshold), threshold = 0; end
    if nargin < 5 || isempty(showing), showing = 15; end   
    if nargin < 6 || isempty(plotting), plotting = 15; end
    
    num_trains = length(spikes);
    if num_trains < 2
        if bitand(showing, 2)
            disp('Not enough trains to calculate the ISI-distance');
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
            existing_ids = [all_figs.Number];
        else
            existing_ids = [];
        end
        while any(existing_ids == fig_profile_id) || any(existing_ids == fig_matrix_id)
                fig_profile_id = fig_profile_id + 2; 
                fig_matrix_id = fig_matrix_id + 2;  
        end
    end 

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
                
                % Avoid null values and division by 0
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

function [MRTS] = autoMRTS(spikes)
% AUTOMRTS Calculates the automated Minimum Resolvable Time Scale (MRTS).
%
%   Computes the Root Mean Square (RMS) of all inter-spike intervals (ISIs) 
%   across all provided spike trains:
%
%   .. math::
%
%      \tau_{MRTS} = \sqrt{ \frac{1}{N_{isi}} \sum_{k=1}^{N_{isi}} (\Delta t_k)^2 }
%
%   where :math:`\Delta t_k` represents the duration of the :math:`k`-th inter-spike interval 
%   and :math:`N_{isi}` is the total number of intervals across all trains.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      MRTS = autoMRTS(spikes);
%
%   :param spikes: Cell array where each element is a vector containing spike timestamps.
%   :type spikes: cell
%
%   :returns: MRTS - Root mean square of inter-spike intervals.
%   :rtype: double
%
%   :Author: Maxime BELTOISE
%   :Date: May 2026


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