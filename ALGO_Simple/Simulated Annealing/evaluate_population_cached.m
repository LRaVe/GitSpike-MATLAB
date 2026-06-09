%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function [P,D,labels,cache] = ...
    evaluate_population_cached( ...
        spikes,...
        neurons,...
        Tmax,...
        Distances,...
        threshold,...
        cache)

    key = mat2str(sort(neurons));

    if isKey(cache,key)

        tmp = cache(key);

        P = tmp.P;
        D = tmp.D;
        labels = tmp.labels;

        return
    end

    [P,D,labels] = evaluate_population( ...
                        spikes,...
                        neurons,...
                        Tmax,...
                        Distances,...
                        threshold);

    tmp.P = P;
    tmp.D = D;
    tmp.labels = labels;

    cache(key) = tmp;

end