%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function spikes = generate_SP_dataset(params)
    
    N = params.N;
    c = params.c;
    
    S = params.S;
    R = params.R;
    
    Tmax = params.Tmax;
    rate = params.rate;
    
    spikes = cell(N,S,R);
    
    %% =========================
    %% STIMULI
    %% =========================
    
    for s = 1:S
    
        %% pooled train du stimulus
    
        nSpikes = poissrnd(rate*Tmax);
    
        pooledTrain = sort(rand(1,nSpikes)*Tmax);
    
        %% repetitions
    
        for r = 1:R
    
            %% reset coding neurons
    
            for n = 1:c
                spikes{n,s,r} = [];
            end
    
            %% -------------------------
            %% coding neurons
            %% -------------------------
    
            for k = 1:length(pooledTrain)
    
                neuron = randi(c);
    
                spikes{neuron,s,r} = [spikes{neuron,s,r} pooledTrain(k)];
    
            end
    
            %% -------------------------
            %% non-coding neurons
            %% -------------------------
    
            for n = c+1:N
    
                nNoise = poissrnd(rate*Tmax);
    
                spikes{n,s,r} = sort(rand(1,nNoise)*Tmax);
    
            end
        end
    end
end


