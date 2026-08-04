%% ALGO Computation
% Author: Maxime BELTOISE
% Date: July 2026


function plot_LL_results(spikes,result)
% PLOT_LL_RESULTS Displays diagnostic figure windows for Labeled Line analysis.
%
%   Generates independent figure windows for detailed visual inspection of 
%   spike train matrices, single-unit distance matrices, discrimination matrices, 
%   performance matrices, and overall population metrics.
%
%   Valid call structure:
%
%   .. code-block:: matlab
%
%      plot_LL_results(spikes, result);
%
%   :param spikes: 3D cell array of size `[N x S x R]` containing spike times.
%   :type spikes: cell
%   :param result: Structure produced by :func:`evaluate_LL_population`.
%   :type result: struct

    %% =========================================================
    %% STRUCTURED SPIKE MATRIX
    %% =========================================================

    figure('Name','Spike Matrix','Color','w');
    
    [N, S, R] = size(spikes);
    
    tiledlayout(N, S, 'Padding','compact', 'TileSpacing','compact');
    
    maxT = 0;
    
    % first pass to get time max
    for n = 1:N
        for s = 1:S
            for r = 1:R
                if ~isempty(spikes{n,s,r})
                    maxT = max(maxT, max(spikes{n,s,r}));
                end
            end
        end
    end
    
    for n = 1:N
        
        for s = 1:S
            
            nexttile((n-1)*S + s);
            hold on;
            
            for r = 1:R
                
                st = spikes{n,s,r};
                
                if isempty(st)
                    continue;
                end
                
                % each repetition is a small vertical offset
                y = r * ones(size(st));
                
                plot(st, y, 'k|', 'MarkerSize', 4);
            end
            
            % styling per cell
            xlim([0 maxT]);
            ylim([0 R+1]);
            
            if n == 1
                title(sprintf('Stimuli%d', s));
            end

            if n == N
                xlabel(sprintf('Time'));
            end
            
            if s == 1
                ylabel(sprintf('Neuron%d', n));
            end
            
            set(gca,'XTick',[],'YTick',[]);
            box on;
            
            hold off;
        end
    end
    
    sgtitle('Spike trains Matrix');



    %% ================================================================
    %% Pairwise distance matrices D_n (SPIKE distance)
    %% ================================================================
    
    figure('Name','Pairwise distance matrices','Color','w');
    
    for n = 1:N

        subplot(2, ceil(N/2), n);
        
        imagesc(result.DistanceMatrix{n});
        axis square;
        colormap([0,0,0;jet]);
        colorbar;
        
        title(sprintf('Neuron %d', n));
        xlabel('Trials');
        ylabel('Trials');
    end
    
    sgtitle('Pairwise SPIKE distance matrices D_n');



    %% ============================================================
    %% Discrimination matrices
    %% ============================================================
    
    Colormap = [0,0,0;hsv(N)];
    
    figure('Name','Discrimination matrices','Color','w');
    
    for n=1:N
    
        subplot(2,ceil(N/2),n)
    
        imagesc(result.Discrimination{n}*n)
    
        axis square
    
        colormap(gca,Colormap)
    
        caxis([0 N])
    
        xticks(1:S)
        yticks(1:S)
    
        xlabel('Stimulus')
        ylabel('Stimulus')
    
        title(sprintf('M_{%d}',n))
    
        DrawSeparators(S,1,1.5);
    
    end
    
    sgtitle('Discrimination matrices')



    %% ============================================================
    %% PANEL E - Performance matrices
    %% ============================================================
    
    figure('Name','Performance matrices');
    
    maxPerformance = 0;
    
    for n = 1:N
        maxPerformance = max(maxPerformance, max(result.Mn{n}(:)));
    end
    
    for n = 1:N
    
        subplot(2,ceil(N/2),n)
    
        imagesc(result.Mn{n})
        axis square
    
        colormap([0,0,0;jet])
    
        caxis([0 maxPerformance])
    
        colorbar
    
        title(sprintf('P_{%d}',n))
    
        DrawSeparators(S,1,1.5);
    
    end
    
    sgtitle('Performance matrices')



    %% ============================================================
    %% Population performance matrix
    %% ============================================================
    
    figure('Name','Population Performance','Color','w');
    
    imagesc(result.populationPerformance);
    
    axis square;
    
    colormap([0,0,0;jet]);
    
    colorbar;
    
    title(sprintf('Population performance (PLL = %.3f)',result.bestP),'FontWeight','bold');
    
    xlabel('Stimulus');
    ylabel('Stimulus');
    
    xticks(1:S);
    yticks(1:S);
    
    DrawSeparators(S,1,1.5);


    %% ============================================================
    %% Best neuron matrix
    %% ============================================================
    
    Colormap = [0,0,0;hsv(N)];

    figure('Name','Best neuron','Color','w');
    
    imagesc(result.bestNeuronMatrix);
    
    axis square;
    
    colormap(gca,Colormap);
    
    caxis([0 N]);
    
    xticks(1:S);
    yticks(1:S);
    
    xlabel('Stimulus');
    ylabel('Stimulus');
    
    title(sprintf('Best neurons  [%s]',sprintf('%d ',result.bestPopulation)),'FontWeight','bold');
    
    DrawSeparators(S,1,1.5);
    
   
    %% ===========================================================
    %% Print summary
    %% ===========================================================
    
    fprintf('\n');
    fprintf('=========== LL RESULTS ===========\n');
    fprintf('Best population : ');
    fprintf('%d ',result.bestPopulation);
    fprintf('\n');
    fprintf('Performance PLL : %.4f\n',result.bestP);
    fprintf('==================================\n');

end