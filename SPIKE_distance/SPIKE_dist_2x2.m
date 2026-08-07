%% ALGO Computation
% Author: Maxime BELTOISE
% Date: July 2026
% Function use to choose the mex or matlab version

function [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold)
% SPIKE_DIST_2X2 Bivariate SPIKE-distance dispatcher.
%
%   Routes the bivariate SPIKE-distance calculation between two spike trains 
%   to either a compiled C++/MEX implementation (`SPIKE_dist_2x2_mex`) or a pure 
%   MATLAB implementation (`SPIKE_dist_2x2_matlab`) based on the global `useMex` flag.
%
%   Syntax:
%      [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2(spikes1, spikes2, ...
%          t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold)
%
%   :param spikes1: Spike timestamps for the first train.
%   :type spikes1: double
%   :param spikes2: Spike timestamps for the second train.
%   :type spikes2: double
%   :param t_min: Start boundary of the observation window.
%   :type t_min: double
%   :param t_max: End boundary of the observation window.
%   :type t_max: double
%   :param aux1_begin: Flag indicating if first spike of train 1 is auxiliary.
%   :type aux1_begin: logical | integer
%   :param aux1_end: Flag indicating if last spike of train 1 is auxiliary.
%   :type aux1_end: logical | integer
%   :param aux2_begin: Flag indicating if first spike of train 2 is auxiliary.
%   :type aux2_begin: logical | integer
%   :param aux2_end: Flag indicating if last spike of train 2 is auxiliary.
%   :type aux2_end: logical | integer
%   :param Distances: Logical 1x4 mask selecting measures: [SPIKE, RI-SPIKE, A-SPIKE, RIA-SPIKE].
%   :type Distances: double | logical
%   :param threshold: Adaptive threshold value for A-SPIKE and RIA-SPIKE variants.
%   :type threshold: double
%
%   :returns: 
%       * **SPIKE_distance_2x2** (*1x4 double*) -- Overall distance values for requested measures.
%       * **profile_mat** (*1x4 cell*) -- Discontinuous time profiles `[t, distance]` for each measure.
%
%   .. important::
%      Routes execution to `SPIKE_dist_2x2_mex` when the global variable `useMex` is `true`,
%      otherwise defaults to `SPIKE_dist_2x2_matlab`.

    global useMex

    if useMex
        [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2_mex(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold);
    else
        [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2_matlab(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold);
    end
end