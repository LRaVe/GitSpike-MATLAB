%% ISI-distance adaptive computationwith auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026

function [dist_matrix,I, I_mean] = f_ISI_distance_adaptive(spikes_trains, tmin, tmax, MRTS)
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
    
    % Manages the MRTS parameter
    if nargin < 4 || isempty(MRTS)
        MRTS = 0; 
    end
    
    if ischar(MRTS) && strcmpi(MRTS, 'auto')
        MRTS = autoMRTS(spikes_trains,'auto'); 
        mode_label = sprintf('Adaptive (auto MRTS = %.3f)', MRTS);
    elseif MRTS > 0
        mode_label = sprintf('Adaptive (manual MRTS = %.3f)', MRTS);
    else
        mode_label = 'Classic (MRTS = 0)';
    end

    num_trains = length(spikes_trains);
    I = [];
    dist_matrix = zeros(num_trains, num_trains);
    all_t_events = [tmin, tmax];
    pair_data = {}; 
    
    % Edge correction
    spikes = cell(1,num_trains);
    for i = 1:num_trains
        s = unique(spikes_trains{i}); 
        spikes{i} = s(s > tmin & s < tmax); 
    end

    if num_trains >= 2
   compteur = 1;
   num_pairs = (num_trains * (num_trains - 1)) / 2;
   num_cols = 2; % On fixe 2 colonnes pour que ce soit lisible
   num_rows = ceil(num_pairs / num_cols);
   figure('Name', ['ISI Evolution - ' mode_label]);
   for i = 1:num_trains
       for j = i+1:num_trains
           compteur = compteur + 1;
           t_all = unique([tmin, spikes{i}, spikes{j}, tmax]);
           all_t_events = [all_t_events, t_all];
           Iij = 0;
           It_list = [];
           for k = 1 : length(t_all)-1
               t_mid = (t_all(k) + t_all(k+1)) / 2;

               %% Correction edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % For the train i
               if isempty(spikes{i})
                   val_x = tmax - tmin; 
               elseif t_mid < spikes{i}(1)
                   % Equation (A.1) 
                   if length(spikes{i}) >= 2
                       val_x = max(spikes{i}(1) - tmin, spikes{i}(2) - spikes{i}(1));
                   else
                       val_x = spikes{i}(1) - tmin;
                   end
               elseif t_mid > spikes{i}(end)
                   % Equation (A.2)
                   if length(spikes{i}) >= 2
                       val_x = max(tmax - spikes{i}(end), spikes{i}(end) - spikes{i}(end-1));
                   else
                       val_x = tmax - spikes{i}(end);
                   end
               else
                   % Intervalle standard (Équation 3)
                   idx = find(spikes{i} <= t_mid, 1, 'last');
                   val_x = spikes{i}(idx+1) - spikes{i}(idx);
               end

               % For the train j
               if isempty(spikes{j})
                   val_y = tmax - tmin;
               elseif t_mid < spikes{j}(1)
                   if length(spikes{j}) >= 2
                       val_y = max(spikes{j}(1) - tmin, spikes{j}(2) - spikes{j}(1));
                   else
                       val_y = spikes{j}(1) - tmin;
                   end
               elseif t_mid > spikes{j}(end)
                   if length(spikes{j}) >= 2
                       val_y = max(tmax - spikes{j}(end), spikes{j}(end) - spikes{j}(end-1));
                   else
                       val_y = tmax - spikes{j}(end);
                   end
               else
                   idy = find(spikes{j} <= t_mid, 1, 'last');
                   val_y = spikes{j}(idy+1) - spikes{j}(idy);
               end                

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               I_t = abs(val_x - val_y) / max([val_x, val_y,MRTS]);
               Iij = Iij+I_t * (t_all(k+1) - t_all(k));
               It_list = [It_list, I_t];
           end
           dist_matrix(i,j) = Iij / (tmax - tmin);
           dist_matrix(j,i) = Iij / (tmax - tmin);
           I = [I,(1/(tmax-tmin))*Iij];
           % display (It_list);
           % display (t_all);
           I_plot = [It_list, It_list(end)];

           pair_data{end+1}.t = t_all;
           pair_data{end}.It = It_list;

           subplot(num_rows, num_cols, compteur - 1); 
           stairs(t_all, I_plot); 
           xlabel('Time');
           ylabel('I_t');
           xlim([0 tmax]);   
           ylim([0 1]);
           title(['Pair ', num2str(i), ' & ', num2str(j)]);  
           subtitle({['Distance: ', num2str(Iij/(tmax-tmin), '%.3f')], ...
               ['Mode: ' mode_label]});
           grid on;
       end
   end   
       I_mean = mean(I);
       t_global = unique(all_t_events); 
       I_matrix = zeros(length(pair_data), length(t_global)-1);
       
       for p = 1:length(pair_data)
           t_p = pair_data{p}.t;
           It_p = pair_data{p}.It;
           for k = 1:length(t_global)-1
               t_mid = (t_global(k) + t_global(k+1)) / 2;
               idx = find(t_p(1:end-1) <= t_mid, 1, 'last');
               I_matrix(p, k) = It_p(idx);
           end
       end
       I_pop_mean = mean(I_matrix, 1);
       
       %display(I)
       
       figure('Name', ['Matrix - ' mode_label]);
       imagesc(dist_matrix); 
       colorbar;
       title(['ISI-distance Matrix - ' mode_label]);
       
       figure('Name', ['Population Average - ' mode_label]);
       stairs(t_global, [I_pop_mean, I_pop_mean(end)]);
       xlabel('Time'); 
       ylabel('Average I_t');
       xlim([0 tmax]);   
       ylim([0 1]);
       title(['Population Average ISI distance - ' mode_label]);
    else
       I_mean = 0; % Cas avec moins de 2 trains
    end
end


function [MRTS] = autoMRTS(spikes, threshold)
    if ischar(threshold) && strcmpi(threshold, 'auto')
        sum_isi_sqr = 0;
        num_isi = 0;
        for i=1:length(spikes)
            for j=1:(length(spikes{i})-1)
                sum_isi_sqr = sum_isi_sqr + (spikes{i}(j+1)-spikes{i}(j))^2;
                num_isi = num_isi + 1;
            end
        end
        MRTS = (sum_isi_sqr/num_isi)^0.5;
    else
        MRTS = threshold;
    end
end