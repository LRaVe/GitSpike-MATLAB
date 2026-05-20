 %% ISI-distance adaptive computationwith auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026


function [dist_matrix,I, I_mean] = f_ISI_distance_adaptive_v1(spikes_trains, ...
    tmin, tmax, threshold,showing,plotting)
% COMPUTE_ADAPTIVE_ISI Calculates the Adaptive ISI-distance between spike trains
%
% Inputs:
%   spikes_trains : cell array containing spike timing vectors
%   tmin, tmax    : time boundaries (minimum and maximum)
%   MRST
%
% Outputs:
%   dist_matrix   : matrix of pairwise adaptive ISI-distances
%   A_I_mean      : average population adaptive ISI-distance
%   (Generates plots showing the pairwise adaptive ISI distance)

       % 1. Manages of the MRTS parameter
    if nargin < 4 || isempty(threshold)
        threshold = 0; 
    end
    if nargin < 5 || isempty(showing)
        showing = 15; 
    end   % Tout afficher par défaut (2+4+8)
    if nargin < 6 || isempty(plotting)
        plotting = 15; 
    end

    if isnumeric(threshold) && threshold <= 0
        error('ISI_Distance:ClassicModeNotAllowed', ...
              'Error: Classic ISI-distance (MRTS = 0) is disabled. Please provide a positive threshold or use ''auto''.');
    end

    MRTS = autoMRTS(spikes_trains, threshold);
    mode_label = sprintf('Adaptive (MRTS = %.3f)', MRTS);

    num_trains = length(spikes_trains);
    I = [];
    dist_matrix = zeros(num_trains, num_trains);
    I_mean = 0;
    all_t_events = [tmin, tmax];
    pair_data = {}; 
    
    spikes = cell(1,num_trains);
    for i = 1:num_trains
        s = unique(spikes_trains{i}); 
        spikes{i} = s(s > tmin & s < tmax); 
    end

    if num_trains >= 2
        compteur = 0;
        num_cols = 2; 
        num_rows = ceil((num_trains * (num_trains - 1) / 2) / num_cols);

        if bitand(plotting, 4)
            figure('Name', ['ISI Evolution - ' mode_label]);
        end

        for i = 1:num_trains
            for j = i+1:num_trains
                compteur = compteur + 1;
                t_all = unique([tmin, spikes{i}, spikes{j}, tmax]);
                all_t_events = [all_t_events, t_all];
                It_list = zeros(1, length(t_all)-1);

                for k = 1 : length(t_all)-1
                    t_mid = (t_all(k) + t_all(k+1)) / 2;
                    
                    % Train i
                    if isempty(spikes{i})
                        vx = tmax - tmin;
                    elseif t_mid < spikes{i}(1)
                        vx = spikes{i}(1) - tmin;
                    elseif t_mid > spikes{i}(end)
                        vx = tmax - spikes{i}(end);
                    else
                        idx = find(spikes{i} <= t_mid, 1, 'last');
                        vx = spikes{i}(idx+1) - spikes{i}(idx);
                    end

                    % Train j
                    if isempty(spikes{j})
                        vy = tmax - tmin;
                    elseif t_mid < spikes{j}(1)
                        vy = spikes{j}(1) - tmin;
                    elseif t_mid > spikes{j}(end)
                        vy = tmax - spikes{j}(end);
                    else
                        idy = find(spikes{j} <= t_mid, 1, 'last');
                        vy = spikes{j}(idy+1) - spikes{j}(idy);
                    end

                    % A_ISI_distance
                    It_list(k) = abs(vx - vy) / max([vx, vy, MRTS]);
                end

                Iij = sum(It_list .* diff(t_all)) / (tmax - tmin);
                dist_matrix(i,j) = Iij;
                dist_matrix(j,i) = Iij;
                I = [I, Iij];
                
                pair_data{compteur}.t = t_all;
                pair_data{compteur}.It = It_list;
                
                if bitand(plotting, 4)
                    subplot(num_rows, num_cols, compteur); 
                    stairs(t_all, [It_list, It_list(end)]); 
                    title(['Pair', num2str(i), ' & ', num2str(j)]);  
                    subtitle(['Dist: ', num2str(Iij, '%.4f')]);
                    ylim([0 1]); 
                    grid on;
                end
            end
        end   
        
        I_mean = mean(I);

        if bitand(showing, 2)
            fprintf('The global adaptative ISI-distance is: %.4f\n', I_mean);
        end
        
        % Average of the population
        t_global = unique(all_t_events); 
        I_matrix = zeros(length(pair_data), length(t_global)-1);
        for p = 1:length(pair_data)
            for k = 1:length(t_global)-1
                t_mid = (t_global(k) + t_global(k+1)) / 2;
                idx = find(pair_data{p}.t(1:end-1) <= t_mid, 1, 'last');
                I_matrix(p, k) = pair_data{p}.It(idx);
            end
        end
        I_pop_mean = mean(I_matrix, 1);

        % global plots

        if bitand(plotting, 8)
            figure('Name', ['Matrix - ' mode_label]);
            imagesc(dist_matrix); 
            colorbar; 
            title('ISI Matrix');
            subtitle(['Global adaptative ISI-distance: ', num2str(I_mean, '%.4f')]);
        end

        if bitand(showing, 8)
            disp('Final adaptative ISI-Distance matrix:');
            disp(dist_matrix);
        end

        if bitand(plotting, 4)
            figure('Name', ['Population Average - ' mode_label]);
            stairs(t_global, [I_pop_mean, I_pop_mean(end)]);
            title('Population Average');
            subtitle(['Global adaptative ISI-distance: ', num2str(I_mean, '%.4f')]);
            ylim([0 1]); 
            grid on;
        end

        if bitand(showing, 4)
            fprintf('\n=== Final adaptative ISI-Distance plot (%s) : ===\n', mode_label);
            fprintf('  Time(t)  |  Average ISI Distance I(t)\n');
            I_pop_extended = [I_pop_mean, I_pop_mean(end)];
            for idx_plot = 1:length(t_global)
                fprintf('      %8.4f     |      %8.4f\n', t_global(idx_plot), ...
                    I_pop_extended(idx_plot));
            end
        end
    end
end 

% subfunction for the MRTS parameter
function [MRTS] = autoMRTS(spikes, threshold)
    if ischar(threshold) && strcmpi(threshold, 'auto')
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
    else
        MRTS = threshold;
    end
    display(MRTS)
end