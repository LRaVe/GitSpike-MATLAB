%% ALGO Computation - Mex version
% Author: Maxime BELTOISE
% Date: July 2026


function P = compute_discrimination_performance_mex(CellMatrix, mask, num_stimuli, num_repetitions, t1, t2, Distances)
    neuron_ids = find(mask == 1);
    
    if isempty(neuron_ids)
        P = 0;
        return;
    end

    [trials, labels] = build_trials(CellMatrix, neuron_ids);

    Tmax = t2; 
    D = compute_population_distance_matrix(trials, Tmax, Distances, 'auto');
    P = mex_compute_discrimination_performance(D, labels);
end