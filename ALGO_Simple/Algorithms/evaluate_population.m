%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function [P,D,labels] = evaluate_population(spikes,neurons,Tmax,Distances,threshold)

    %% pooled trials

    [trials,labels] = build_trials(spikes,neurons);

    %% distance matrix

    D = compute_population_distance_matrix(trials,Tmax,Distances,threshold);

    %% performance

    P = compute_discrimination_performance(D,labels);

end