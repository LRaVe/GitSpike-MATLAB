%% ISI-distance adaptive computationwith auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026

function [dist_matrix,I, I_mean] = f_ISI_distance_adaptive_v1(spikes_trains, tmin, tmax, threshold)
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
        MRTS = 0; 
    else
        MRTS = autoMRTS(spikes_trains, threshold); 
    end
    
    if MRTS > 0
        mode_label = sprintf('Adaptive (MRTS = %.3f)', MRTS);
    else
        mode_label = 'Classic (MRTS = 0)';
    end

    num_trains = length(spikes_trains);
    I = [];
    dist_matrix = zeros(num_trains, num_trains);
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
        figure('Name', ['ISI Evolution - ' mode_label]);

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

                subplot(num_rows, num_cols, compteur); 
                stairs(t_all, [It_list, It_list(end)]); 
                title(['Pair ', num2str(i), ' & ', num2str(j)]);  
                subtitle(['Dist: ', num2str(Iij, '%.4f')]);
                ylim([0 1]); grid on;
            end
        end   
        
        I_mean = mean(I);
        
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
        figure('Name', ['Matrix - ' mode_label]);
        imagesc(dist_matrix); colorbar; title('ISI Matrix');
        
        figure('Name', ['Population Average - ' mode_label]);
        stairs(t_global, [I_pop_mean, I_pop_mean(end)]);
        title('Population Average'); ylim([0 1]); grid on;

    else
        I_mean = 0;
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