%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function D = compute_population_distance_matrix(trials,Tmax)
    
    T = length(trials);
    
    D = zeros(T);
    
    %% =====================================
    %% auxiliary spikes
    %% =====================================
    
    [trials_aux,aux_begin,aux_end] = ...
        add_auxiliary_spikes(trials,0,Tmax);
    
    %% =====================================
    %% pairwise distances
    %% =====================================
    
    Distances = [1 0 0 0];
    threshold = 0;
    
    for i = 1:T
    
        for j = i+1:T
    
            [d,~] = SPIKE_dist_2x2( ...
                trials_aux{i}, ...
                trials_aux{j}, ...
                0, Tmax, ...
                aux_begin(i), ...
                aux_end(i), ...
                aux_begin(j), ...
                aux_end(j), ...
                Distances, ...
                threshold);
    
            D(i,j) = d(1);
            D(j,i) = d(1);
    
        end
    end
end


