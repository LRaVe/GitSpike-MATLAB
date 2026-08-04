%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026



%% =========================================================
% Windowing and edge correction
% =========================================================


function [spikes, aux_begin, aux_end] = add_auxiliary_spikes(spikes, t_min, t_max)
    % Window spike trains and add auxiliary spike times at the
    % beginning and end of the observation window when needed.

    % ===============================================
    % ====== Handle cell array of spike trains ======
    % ===============================================

    if iscell(spikes)
        aux_begin = zeros(1, length(spikes));
        aux_end   = zeros(1, length(spikes));
        for i = 1:length(spikes)
            [spikes{i}, aux_begin(i), aux_end(i)] = process_single_train(spikes{i}, t_min, t_max);
        end
        return;
    end

    % =======================================
    % ====== Handle single spike train ======
    % =======================================

    [spikes, aux_begin, aux_end] = process_single_train(spikes, t_min, t_max);

end


function [train, aux_begin, aux_end] = process_single_train(train, t_min, t_max)

    aux_begin = 0;
    aux_end   = 0;

    train = sort(unique(train));

    % =======================================
    % Empty train
    % =======================================

    if isempty(train)
        train = [t_min t_max];
        aux_begin = 1;
        aux_end   = 1;
        return;
    end

    % =======================================
    % Windowing
    % =======================================

    idx_before = find(train < t_min, 1, 'last');
    idx_after  = find(train > t_max, 1, 'first');

    idx_inside = find(train >= t_min & train <= t_max);

    new_train = [];

    % Keep last spike before t_min
    if ~isempty(idx_before)
        new_train(end+1) = train(idx_before);
    end

    % Keep spikes inside window
    if ~isempty(idx_inside)
        new_train = [new_train train(idx_inside)];
    end

    % Keep first spike after t_max
    if ~isempty(idx_after)
        new_train(end+1) = train(idx_after);
    end

    train = new_train;

    % =======================================
    % Left edge correction
    % =======================================

    if train(1) > t_min

        if length(train) >= 2
            aux = train(1) - max(train(1)-t_min, train(2)-train(1));
        else
            aux = t_min;
        end

        train = [aux train];
        aux_begin = 1;
    end

    % =======================================
    % Right edge correction
    % =======================================

    if train(end) < t_max

        if length(train) >= 2
            aux = train(end) + max(t_max-train(end), train(end)-train(end-1));
        else
            aux = t_max;
        end

        train = [train aux];
        aux_end = 1;
    end

end

