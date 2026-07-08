%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function [P,D,labels,cache] = evaluate_population_cached(spikes,neurons,Tmax,Distances,threshold,cache)

    %% =====================================
    %% unique key of the population
    %% =====================================

    key = mat2str(sort(neurons));

    %% =====================================
    %% cache hit
    %% =====================================

    if isKey(cache.data,key)

        cache.hits = cache.hits + 1;

        tmp = cache.data(key);

        P = tmp.P;
        D = tmp.D;
        labels = tmp.labels;

        return

    end

    %% =====================================
    %% cache miss
    %% =====================================

    cache.misses = cache.misses + 1;

    [P,D,labels] = evaluate_population(spikes,neurons,Tmax,Distances,threshold);

    tmp.P = P;
    tmp.D = D;
    tmp.labels = labels;

    cache.data(key) = tmp;

end