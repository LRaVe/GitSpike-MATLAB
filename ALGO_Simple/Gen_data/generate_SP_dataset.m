%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function spikes = generate_SP_dataset(params)
% GENERATE_SP_DATASET Generates synthetic spike trains for Summed Population analysis.
%
%   Generates a 3D cell array of simulated spike trains based on input 
%   neural population parameters. The generated dataset models multi-neuron 
%   responses across multiple stimuli and trial repetitions for Summed Population (SP) 
%   benchmarks, including Collective (COLL), Individual (INDI), and Non-Coding neurons.
%
%   Valid call structure:
%
%   .. code-block:: matlab
%
%      spikes = generate_SP_dataset(params);
%
%   :param params: Structure containing dataset generation parameters:
%                  
%                  * **N** (*integer*): Total number of neurons in the population.
%                  * **c** (*integer*): Number of collective (COLL) neurons.
%                  * **nIndi** (*integer*): Number of individual (INDI) neurons.
%                  * **S** (*integer*): Number of distinct stimulus conditions.
%                  * **R** (*integer*): Number of trial repetitions per stimulus.
%                  * **Tmax** (*double*): Upper temporal boundary of trial duration (in seconds).
%                  * **rate** (*double*): Baseline firing rate (in Hz).
%                  * **indiJitter** (*double*): Gaussian jitter standard deviation for INDI neurons.
%   :type params: struct
%
%   :returns: 
%             * **spikes** (*cell*) -- 3D cell array of dimensions `[N x S x R]` containing spike timestamps.

    N    = params.N;

    c    = params.c;      % neurones Coll
    nIndi = params.nIndi; % neurones Indi

    S    = params.S;
    R    = params.R;

    Tmax = params.Tmax;
    rate = params.rate;
    indiJitter = params.indiJitter;

    spikes = cell(N,S,R);

    %% =====================================================
    %% COLL neurons
    %% =====================================================

    pooledRate = c*rate;

    for s = 1:S

        nSpikes = poissrnd(pooledRate*Tmax);

        pooledTrain = sort(rand(1,nSpikes)*Tmax);

        for r = 1:R

            assignment = randi(c,1,nSpikes);

            for n = 1:c

                spikes{n,s,r} = pooledTrain(assignment==n);

            end
        end
    end

    %% =====================================================
    %% INDI neurons
    %% =====================================================

    for n = c+1:c+nIndi

        template = cell(1,S);
    
        for s = 1:S
    
            nSpikes = poissrnd(rate*Tmax);
    
            template{s} = sort(rand(1,nSpikes)*Tmax);
    
        end
    
        for s = 1:S
    
            for r = 1:R
    
                train = template{s};
    
                train = train + indiJitter*randn(size(train));
    
                train(train<0) = 0;
                train(train>Tmax) = Tmax;
    
                spikes{n,s,r} = sort(train);
    
            end
    
        end
    
    end

    %% =====================================================
    %% NON CODING
    %% =====================================================

    for n = c+nIndi+1:N

        for s = 1:S

            for r = 1:R

                nNoise = poissrnd(rate*Tmax);

                spikes{n,s,r} = sort(rand(1,nNoise)*Tmax);

            end

        end

    end

end



