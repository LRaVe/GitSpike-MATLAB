%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026

function [D_global, profile_global, D_matrix, RI_D_global, RI_profile_global, RI_D_matrix, A_D_global, A_profile_global, A_D_matrix, RIA_D_global, RIA_profile_global, RIA_D_matrix] = SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end, Distances, threshold)
% SPIKE_DIST_N Computes multivariate SPIKE-distance across N spike trains.
%
%   Evaluates all pairwise combinations of N spike trains, builds the pairwise 
%   distance matrices, and constructs global ensemble profiles and scalar averages.
%
%   Syntax:
%      [D_global, profile_global, D_matrix, ...
%       RI_D_global, RI_profile_global, RI_D_matrix, ...
%       A_D_global, A_profile_global, A_D_matrix, ...
%       RIA_D_global, RIA_profile_global, RIA_D_matrix] = ...
%          SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end, Distances, threshold)
%
%   :param spikes: Cell array containing spike timestamp vectors for each train.
%   :type spikes: cell
%   :param t_min: Lower bound of the evaluation window.
%   :type t_min: double
%   :param t_max: Upper bound of the evaluation window.
%   :type t_max: double
%   :param aux_begin: Boolean array indicating auxiliary start spikes per train.
%   :type aux_begin: logical | double
%   :param aux_end: Boolean array indicating auxiliary end spikes per train.
%   :type aux_end: logical | double
%   :param Distances: 1x4 mask selecting enabled distance variants `[SPIKE, RI, A, RIA]`.
%   :type Distances: double | logical
%   :param threshold: Adaptive distance calculation threshold.
%   :type threshold: double
%
%   :returns:
%       * **D_global** (*double*) -- Global mean SPIKE-distance across all pairs.
%       * **profile_global** (*Nx2 double*) -- Averaged time profile for SPIKE-distance.
%       * **D_matrix** (*NxN double*) -- Pairwise distance matrix for SPIKE-distance.
%       * **RI_D_global, RI_profile_global, RI_D_matrix** -- Corresponding outputs for RI-SPIKE-distance.
%       * **A_D_global, A_profile_global, A_D_matrix** -- Corresponding outputs for A-SPIKE-distance.
%       * **RIA_D_global, RIA_profile_global, RIA_D_matrix** -- Corresponding outputs for RIA-SPIKE-distance.

    N = length(spikes);

    D_matrix = zeros(N);
    RI_D_matrix = zeros(N);
    A_D_matrix = zeros(N);
    RIA_D_matrix = zeros(N);

    profiles = cell(1,4);
    idx_prof = 0;

    %% =====================================================
    % PAIRWISE DISTANCES
    %% =====================================================

    for i = 1:N
        for j = i+1:N

            [d, prof] = SPIKE_dist_2x2(spikes{i}, spikes{j}, t_min, t_max, aux_begin(i), aux_end(i), aux_begin(j), aux_end(j), Distances, threshold);
            if numel(d) >= 1
                D_matrix(i,j) = d(1);
                D_matrix(j,i) = d(1);
            end
            if numel(d) >= 2
                RI_D_matrix(i,j) = d(2);
                RI_D_matrix(j,i) = d(2);
            end
            if numel(d) >= 3
                A_D_matrix(i,j) = d(3);
                A_D_matrix(j,i) = d(3);
            end
            if numel(d) >= 4
                RIA_D_matrix(i,j) = d(4);
                RIA_D_matrix(j,i) = d(4);
            end

            idx_prof = idx_prof + 1;

            profiles{1}{idx_prof} = prof{1};    % Classic
            profiles{2}{idx_prof} = prof{2};    % RI
            profiles{3}{idx_prof} = prof{3};    % A
            profiles{4}{idx_prof} = prof{4};    % RIA
        end
    end

    %% =====================================================
    % GLOBAL DISTANCE
    %% =====================================================

    D_global = mean(D_matrix(triu(true(N),1)));
    RI_D_global = mean(RI_D_matrix(triu(true(N),1)));
    A_D_global = mean(A_D_matrix(triu(true(N),1)));
    RIA_D_global = mean(RIA_D_matrix(triu(true(N),1)));

    %% =====================================================
    % ALL TIME COORDINATES
    %% =====================================================

    %N = number_spikes*2*(factorial(num_trains)/(2*factorial(num_trains-2)));
    
    base_profile_idx = find(Distances, 1, 'first');
    t_all = zeros(sum(cellfun(@length,profiles{base_profile_idx})),1);

    idx = 1;

    if ~isempty(base_profile_idx)
        for p = 1:length(profiles{base_profile_idx})
            if ~isempty(profiles{base_profile_idx}{p})
                v = profiles{base_profile_idx}{p}(:,1);
                n = length(v);
    
                t_all(idx:idx+n-1) = v;
                idx = idx + n;
            end
        end
    end

    t_all = unique(sort(t_all));

    %% =====================================================
    % GLOBAL PROFILE SPIKE-distance
    %% =====================================================

    prof = cell(1,4);
    for i=1:4
        
        for k = 1:length(t_all)
    
            t = t_all(k);
    
            vals_left  = zeros(1,length(profiles{i}));
            vals_right = zeros(1,length(profiles{i}));
    
            has_discontinuity = false;
    
            %% -------------------------------------------------
            % scan all pairwise profiles
            %% -------------------------------------------------
    
            for p = 1:length(profiles{i})
    
                P = profiles{i}{p};

                if isempty(P)
                    continue;
                end
    
                idx = find(P(:,1) == t);
    
                %% =============================================
                % CASE 1 : discontinuity in this profile
                %% =============================================
    
                if length(idx) == 2
    
                    has_discontinuity = true;
    
                    vals_left(p)  = P(idx(1),2);
                    vals_right(p) = P(idx(2),2);
    
                %% =============================================
                % CASE 2 : single point
                %% =============================================
    
                elseif isscalar(idx) %faster than length(idx) == 1
    
                    vals_left(p)  = P(idx,2);
                    vals_right(p) = P(idx,2);
    
                %% =============================================
                % CASE 3 : interpolation
                %% =============================================
    
                else
    
                    idx_before = find(P(:,1) < t, 1, 'last');
                    idx_after  = find(P(:,1) > t, 1, 'first');
    
                    if ~isempty(idx_before) && ~isempty(idx_after)
    
                        t1 = P(idx_before,1);
                        t2 = P(idx_after,1);
    
                        S1 = P(idx_before,2);
                        S2 = P(idx_after,2);
    
                        %% avoid division by zero
                        if t2 ~= t1
                            S_interp = S1 + (S2-S1)*(t-t1)/(t2-t1);
                        else
                            S_interp = S1;
                        end
    
                        vals_left(p)  = S_interp;
                        vals_right(p) = S_interp;
    
                    end
                end
            end

            %% -------------------------------------------------
            % averaging
            %% -------------------------------------------------
    
            if has_discontinuity
                prof{i} = [ prof{i} ; t mean(vals_left) ; t mean(vals_right)];
            else
                prof{i} = [prof{i} ; t mean(vals_left)];
            end
        end
    end
    
    % sort
    profile_global = sortrows(prof{1},1);
    RI_profile_global = sortrows(prof{2},1);
    A_profile_global = sortrows(prof{3},1);
    RIA_profile_global = sortrows(prof{4},1);

end