%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function pooled = pool_neurons(spikes,neurons,s,r)

    %% space allocation

    totalSpikes = 0;

    for n = neurons
        totalSpikes = totalSpikes + numel(spikes{n,s,r});
    end

    pooled = zeros(1,totalSpikes);

    %% filling

    idx = 1;

    for n = neurons

        v = spikes{n,s,r};

        l = numel(v);

        pooled(idx:idx+l-1) = v;

        idx = idx + l;

    end

    pooled = sort(pooled);

end