%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function [trials,labels] = build_trials(spikes,neurons)

    [~,S,R] = size(spikes);

    T = S*R;

    trials = cell(1,T);
    labels = zeros(1,T);

    idx = 0;

    for s = 1:S

        for r = 1:R

            idx = idx + 1;

            trials{idx} = pool_neurons(spikes,neurons,s,r);

            labels(idx) = s;

        end

    end

end

