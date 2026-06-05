%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function spikes = generate_SP_dataset(params)

    N    = params.N;
    c    = params.c;

    S    = params.S;
    R    = params.R;

    Tmax = params.Tmax;
    rate = params.rate;

    spikes = cell(N,S,R);

    pooledRate = c * rate;

    for s = 1:S

        %% train poolé du stimulus

        nSpikes = poissrnd(pooledRate*Tmax);

        pooledTrain = sort(rand(1,nSpikes)*Tmax);

        for r = 1:R

            %% reset neurones codants

            for n = 1:c
                spikes{n,s,r} = [];
            end

            %% distribution aléatoire

            assignment = randi(c,1,nSpikes);

            for n = 1:c

                spikes{n,s,r} = pooledTrain(assignment==n);

            end

            %% neurones non codants

            for n = c+1:N

                nNoise = poissrnd(rate*Tmax);

                spikes{n,s,r} = sort(rand(1,nNoise)*Tmax);

            end

        end

    end

end


