%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function T0 = initialize_temperature( ...
                        spikes,...
                        Tmax,...
                        Distances,...
                        threshold,...
                        N0)

    N = size(spikes,1);

    %% population initiale

    pop = randperm(N,randi([1 N]));

    [Pprev,~,~] = evaluate_population( ...
                        spikes,...
                        pop,...
                        Tmax,...
                        Distances,...
                        threshold);

    deltaP = zeros(1,N0);

    %% marche aléatoire

    for k = 1:N0

        pop = random_neighbor(pop,N);

        [Pnew,~,~] = evaluate_population( ...
                        spikes,...
                        pop,...
                        Tmax,...
                        Distances,...
                        threshold);

        deltaP(k) = abs(Pnew-Pprev);

        Pprev = Pnew;

    end

    meanDelta = mean(deltaP);

    T0 = -meanDelta/log(0.95);

end