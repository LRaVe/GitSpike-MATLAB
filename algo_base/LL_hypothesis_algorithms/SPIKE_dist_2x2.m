%% ALGO Computation
% Author: Maxime BELTOISE
% Date: July 2026
% Function use to choose the mex or matlab version

function [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold)

    global useMex

    if useMex
        [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2_mex(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold);
    else
        [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2_matlab(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold);
    end
end