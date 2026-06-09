%% SPIKE-order pairwise computation
% Author: Agathe JULIEN
% Date: May 2026



function res=pairwise_order(tmin,tmax,spikes,spike_ind1,spike_ind2)
    % Compute pairwise spike ordering between two spike trains
    addpath('spike_common');
    n=length(spikes);
    if spike_ind1>n || spike_ind2>n || spike_ind1<1 || spike_ind2<1
        error('Index out of bounds');
    end

    s1=spikes{spike_ind1};
    s2=spikes{spike_ind2};
    s1 = s1(:)';  % force row vectors
    s2 = s2(:)';
    res=zeros(1,length(s1));

    for i=1:length(s1)
        [min_dist, closest_j] = min(abs(s2 - s1(i)));
        window = coincidence_window(tmin, tmax, spikes, spike_ind1, spike_ind2, i, closest_j);
        % fprintf('i=%d, s1(i)=%.2f, best_j=%d, s2(best_j)=%.2f, dist=%.4f, win=%.4f\n', i, s1(i), closest_j, s2(closest_j), min_dist, window);
        % Determine order: within window = 0, otherwise compare spike times
        if min_dist <= window
            if s2(closest_j) > s1(i)
                res(i) = 1;       % s1 leads → +1
            elseif s2(closest_j) < s1(i)
                res(i) = -1;      % s1 trails → -1
            % else exact tie → stays 0
            end
        % else: no coincident spike found → res(i) stays 0
        end
    end
end


