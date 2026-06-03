%% SPIKE-distance computation with auxiliary boundary spikes
% Author: Maxime BELTOISE
% Date: May 2026


%% =====================================================
% DISPLAY
% =====================================================


function SPIKE_distances_all(spikes, t_min, t_max, threshold, measures, adaptive_measures, showing, plotting, aux_begin, aux_end)

    Distances = [0, 0, 0, 0];

    if mod(measures,4)>1
        Distances(1) = 1;
    end
    if mod(measures,8)>3
        Distances(2) = 1;
    end
    if mod(adaptive_measures,4)>1 && ((threshold > 0) || (Distances(1) == 0))
        Distances(3) = 1;
    end
    if mod(adaptive_measures,8)>3 && ((threshold > 0) || (Distances(2) == 0))
        Distances(4) = 1;
    end
    if any(Distances ~= 0)
        [SPIKE_distance, SPIKE_distance_profile, SPIKE_distance_matrix, RI_SPIKE_distance, RI_SPIKE_distance_profile, RI_SPIKE_distance_matrix, A_SPIKE_distance, A_SPIKE_distance_profile, A_SPIKE_distance_matrix, RIA_SPIKE_distance, RIA_SPIKE_distance_profile, RIA_SPIKE_distance_matrix] = SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end, Distances, threshold);
    end

    if mod(adaptive_measures,4)>1 && (Distances(3) == 0)
        A_SPIKE_distance = SPIKE_distance;
        A_SPIKE_distance_profile = SPIKE_distance_profile;
        A_SPIKE_distance_matrix = SPIKE_distance_matrix;
    end
    if mod(adaptive_measures,8)>3 && (Distances(4) == 0)
        RIA_SPIKE_distance = RI_SPIKE_distance;
        RIA_SPIKE_distance_profile = RI_SPIKE_distance_profile;
        RIA_SPIKE_distance_matrix = RI_SPIKE_distance_matrix;
    end


    if mod(measures,4)>1                                                        % SPIKE-distance
        if mod(showing,4)>1
            fprintf('SPIKE-distance: %.4f\n', SPIKE_distance);
        end
        if mod(showing,8)>3
            fprintf('SPIKE-distance profile: \n');
            disp(SPIKE_distance_profile);
        end
        if mod(showing,16)>7
            fprintf('SPIKE-distance matrix: \n');
            disp(SPIKE_distance_matrix);
        end
        if mod(plotting,8)>3
            titleStr = ['SPIKE-distance = ', num2str(SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            area(SPIKE_distance_profile(:,1), SPIKE_distance_profile(:,2));
            xlabel('Time');
            ylabel('SPIKE-distance');
            title(titleStr);
            xlim([t_min t_max]);
            ylim([0 1]);
            colororder([0.5 0.5 1]);
            grid on;
            box on;
        end
        if mod(plotting,16)>7
            titleStr = ['SPIKE-distance = ', num2str(SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            imagesc(SPIKE_distance_matrix);
            colorbar;
            box on;
            colormap jet;
            set(gca,'XTick',1:length(spikes),'YTick',1:length(spikes));
            xlabel('Spike trains');
            ylabel('Spike trains');
            title(titleStr);
        end
    end
    
    if mod(measures,8)>3                                                        % RI-SPIKE-distance
        if mod(showing,4)>1
            fprintf('RI-SPIKE-distance: %.4f\n', RI_SPIKE_distance);
        end
        if mod(showing,8)>3
            fprintf('RI-SPIKE-distance profile: \n');
            disp(RI_SPIKE_distance_profile);
        end
        if mod(showing,16)>7
            fprintf('RI-SPIKE-distance matrix: \n');
            disp(RI_SPIKE_distance_matrix);
        end
        if mod(plotting,8)>3
            titleStr = ['RI-SPIKE-distance = ', num2str(RI_SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            area(RI_SPIKE_distance_profile(:,1), RI_SPIKE_distance_profile(:,2));
            xlabel('Time');
            ylabel('RI-SPIKE-distance');
            title(titleStr);
            xlim([t_min t_max]);
            ylim([0 1]);
            colororder([0.5 0.5 1]);
            grid on;
            box on;
        end
        if mod(plotting,16)>7
            titleStr = ['RI-SPIKE-distance = ', num2str(RI_SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            imagesc(RI_SPIKE_distance_matrix);
            colorbar;
            box on;
            colormap jet;
            set(gca,'XTick',1:length(spikes),'YTick',1:length(spikes));
            xlabel('Spike trains');
            ylabel('Spike trains');
            title(titleStr);
        end
    end


    if mod(adaptive_measures,4)>1                                               % A-SPIKE-distance
        if mod(showing,4)>1
            fprintf('A-SPIKE-distance: %.4f\n', A_SPIKE_distance);
        end
        if mod(showing,8)>3
            fprintf('A-SPIKE-distance profile: \n');
            disp(A_SPIKE_distance_profile);
        end
        if mod(showing,16)>7
            fprintf('A-SPIKE-distance matrix: \n');
            disp(A_SPIKE_distance_matrix);
        end
        if mod(plotting,8)>3
            titleStr = ['A-SPIKE-distance = ', num2str(A_SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            area(A_SPIKE_distance_profile(:,1), A_SPIKE_distance_profile(:,2));
            xlabel('Time');
            ylabel('A-SPIKE-distance');
            title(titleStr);
            xlim([t_min t_max]);
            ylim([0 1]);
            colororder([0.5 0.5 1]);
            grid on;
            box on;
        end
        if mod(plotting,16)>7
            titleStr = ['A-SPIKE-distance = ', num2str(A_SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            imagesc(A_SPIKE_distance_matrix);
            colorbar;
            box on;
            colormap jet;
            set(gca,'XTick',1:length(spikes),'YTick',1:length(spikes));
            xlabel('Spike trains');
            ylabel('Spike trains');
            title(titleStr);
        end
    end


    if mod(adaptive_measures,8)>3                                               % RIA-SPIKE-distance
        if mod(showing,4)>1
            fprintf('RIA-SPIKE-distance: %.4f\n', RIA_SPIKE_distance);
        end
        if mod(showing,8)>3
            fprintf('RIA-SPIKE-distance profile: \n');
            disp(RIA_SPIKE_distance_profile);
        end
        if mod(showing,16)>7
            fprintf('RIA-SPIKE-distance matrix: \n');
            disp(RIA_SPIKE_distance_matrix);
        end
        if mod(plotting,8)>3
            titleStr = ['RIA-SPIKE-distance = ', num2str(RIA_SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            area(RIA_SPIKE_distance_profile(:,1), RIA_SPIKE_distance_profile(:,2));
            xlabel('Time');
            ylabel('RIA-SPIKE-distance');
            title(titleStr);
            xlim([t_min t_max]);
            ylim([0 1]);
            colororder([0.5 0.5 1]);
            grid on;
            box on;
        end
        if mod(plotting,16)>7
            titleStr = ['RIA-SPIKE-distance = ', num2str(RIA_SPIKE_distance)];
            figure('Name', titleStr, 'NumberTitle', 'off');
            imagesc(RIA_SPIKE_distance_matrix);
            colorbar;
            box on;
            colormap jet;
            set(gca,'XTick',1:length(spikes),'YTick',1:length(spikes));
            xlabel('Spike trains');
            ylabel('Spike trains');
            title(titleStr);
        end
    end
end




%% =========================================================
% FUNCTION
% =========================================================


%% =========================================================
% SPIKE-distance between TWO spike trains (2x2 case)
% =========================================================
function [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end, Distances, threshold)


    %% =====================================================
    % INITIALIZE PROFILE
    %% =====================================================

    SPIKE_distance_profile = cell(1,4);

    %% =====================================================
    % MAIN LOOP
    %% =====================================================

    for idx_1 = 1 : length(spikes1)

        if spikes2(1)>spikes1(idx_1)  %In case first spikes are not synchronised
            idx_2=1; 
        elseif spikes2(end)<=spikes1(idx_1) %In case last spikes are not synchronised
            idx_2=length(spikes2)-1;
        else
            idx_2 = find(spikes2(1:end) <= spikes1(idx_1), 1, 'last');
        end

        %% ---------------------------------------------
        % Train 2 contribution
        %% ---------------------------------------------
        ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);

        % nearest-neighbor distances
        delta_tp_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);

        delta_tf_2 = auxiliary_delta(spikes2(idx_2+1), spikes2, spikes1, idx_2+1, aux2_end);

        xp_2 = spikes1(idx_1) - spikes2(idx_2);
        xf_2 = spikes2(idx_2+1) - spikes1(idx_1);

        S_2 = ((delta_tp_2*xf_2) + (delta_tf_2*xp_2)) / ISI_dist_2;

        if idx_1>1
            ISI_dist_1 = spikes1(idx_1) - spikes1(idx_1-1);

            S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_end);
            
            if Distances(1)==1 %SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
                SPIKE_distance_profile{1}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(2)==1 %RI-SPIKE-distance
                S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));
                SPIKE_distance_profile{2}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(3)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(4)==1 %RIA-SPIKE-distance
                S = (S_1 + S_2) / (2*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{4}{end+1} = [spikes1(idx_1) S];
            end
        end

        if idx_1<length(spikes1)
            ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);

            S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);

            if Distances(1)==1 %SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
                SPIKE_distance_profile{1}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(2)==1 %RI-SPIKE-distance
                S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));
                SPIKE_distance_profile{2}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(3)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(4)==1 %RIA-SPIKE-distance
                S = (S_1 + S_2) / (2*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{4}{end+1} = [spikes1(idx_1) S];
            end
        end
    end

    for idx_2 = 1 : length(spikes2)
    
        if spikes1(1)>spikes2(idx_2)  %In case first spikes are not synchronised
            idx_1=1;
        elseif spikes1(end)<=spikes2(idx_2) %In case last spikes are not synchronised
            idx_1=length(spikes1)-1;
        else
            idx_1 = find(spikes1(1:end) <= spikes2(idx_2), 1, 'last');
        end
        

        %% ---------------------------------------------
        % Train 1 contribution
        %% ---------------------------------------------
        ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);

        % nearest-neighbor distances
        delta_tp_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);

        delta_tf_1 = auxiliary_delta(spikes1(idx_1+1), spikes1, spikes2, idx_1+1, aux1_end);

        xp_1 = spikes2(idx_2) - spikes1(idx_1);
        xf_1 = spikes1(idx_1+1) - spikes2(idx_2);

        S_1 = ((delta_tp_1*xf_1) + (delta_tf_1*xp_1)) / ISI_dist_1;

        if idx_2>1
            ISI_dist_2 = spikes2(idx_2) - spikes2(idx_2-1);

            S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_end);

            if Distances(1)==1 %SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
                SPIKE_distance_profile{1}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(2)==1 %RI-SPIKE-distance
                S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));
                SPIKE_distance_profile{2}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(3)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(4)==1 %RIA-SPIKE-distance
                S = (S_1 + S_2) / (2*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{4}{end+1} = [spikes2(idx_2) S];
            end
        end

        if idx_2<length(spikes2)
            ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);

            S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);

            if Distances(1)==1 %SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
                SPIKE_distance_profile{1}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(2)==1 %RI-SPIKE-distance
                S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));
                SPIKE_distance_profile{2}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(3)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(4)==1 %RIA-SPIKE-distance
                S = (S_1 + S_2) / (2*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{4}{end+1} = [spikes2(idx_2) S];
            end
        end
    end

    profile_mat=cell(1,4);
    
    %% =====================================================
    % CONVERT PROFILE SPIKE-distance
    %% =====================================================
    for k=1:4
        if Distances(k)==1
            profile = cell2mat(SPIKE_distance_profile{k}');
            profile = sortrows(profile, 1);

            % ==========================================
            % keep only points inside interval
            % ==========================================
            for i = 1:size(profile,1)
            
                % ==========================================
                % point before t_min
                % ==========================================
                if profile(i,1) < t_min
            
                    % first index whose abscissa >= t_min
                    idx = find(profile(:,1) >= t_min, 1, 'first');
        
                    % linear interpolation
                    profile(i,2) = profile(i,2) + ((profile(idx,2) - profile(i,2)) / (profile(idx,1) - profile(i,1))) * (t_min - profile(i,1));
            
                    % projection onto boundary
                    profile(i,1) = t_min;
            
                % ==========================================
                % point after t_max
                % ==========================================
                elseif profile(i,1) > t_max
            
                    % last index whose abscissa <= t_max
                    idx = find(profile(:,1) <= t_max, 1, 'last');
            
                    % linear interpolation
                    profile(i,2) = profile(idx,2) + ((profile(i,2) - profile(idx,2)) / (profile(i,1) - profile(idx,1))) * (t_max - profile(idx,1));
            
                    % projection onto boundary
                    profile(i,1) = t_max;

                end
            end
            % remove duplicates
            [~, idx] = unique(profile,'rows','stable');
            profile = profile(idx,:);

            % sort
            profile = sortrows(profile,1);

            profile_mat{k}=profile;
        else
            profile_mat{k} = [];
        end
    end

    

    %% =====================================================
    % FINAL DISTANCE
    %% =====================================================
    
    SPIKE_distance_2x2 = zeros(1,4);

    for k=1:4
        if Distances(k)==1
    
            t = profile_mat{k}(:,1);
            S = profile_mat{k}(:,2);
    
            SPIKE_distance_2x2(k) = trapz(t,S)/(t_max-t_min);
        end
    end
end



%% =========================================================
% AUXILIARY DELTA MANAGEMENT
%
% For auxiliary spikes:
% - beginning auxiliary spike inherits delta from right neighbor
% - ending auxiliary spike inherits delta from left neighbor
%% =========================================================
function delta = auxiliary_delta(spike, own_train, other_train, idx, aux_idx)

    % standard nearest-neighbor distance
    delta_std = min(abs(spike - other_train(:)));
    delta = delta_std;

    % auxiliary at beginning
    if (idx == 1) && aux_idx
        delta = min(abs(own_train(2) - other_train(:)));
    end

    % auxiliary at end
    if (idx == length(own_train)) && aux_idx
        delta = min(abs(own_train(end-1) - other_train(:)));
    end
end



%% =========================================================
% SPIKE-distance between N spike trains
% =========================================================
function [D_global, profile_global, D_matrix, RI_D_global, RI_profile_global, RI_D_matrix, A_D_global, A_profile_global, A_D_matrix, RIA_D_global, RIA_profile_global, RIA_D_matrix] = SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end, Distances, threshold)

    N = length(spikes);

    D_matrix = zeros(N);
    RI_D_matrix = zeros(N);
    A_D_matrix = zeros(N);
    RIA_D_matrix = zeros(N);

    profiles = cell(1,4);
    idx_prof = 0;

    %% =====================================================
    % PAIRWISE DISTANCES
    %% =====================================================

    for i = 1:N
        for j = i+1:N

            [d, prof] = SPIKE_dist_2x2(spikes{i}, spikes{j}, t_min, t_max, aux_begin(i), aux_end(i), aux_begin(j), aux_end(j), Distances, threshold);
            if numel(d) >= 1
                D_matrix(i,j) = d(1);
                D_matrix(j,i) = d(1);
            end
            if numel(d) >= 2
                RI_D_matrix(i,j) = d(2);
                RI_D_matrix(j,i) = d(2);
            end
            if numel(d) >= 3
                A_D_matrix(i,j) = d(3);
                A_D_matrix(j,i) = d(3);
            end
            if numel(d) >= 4
                RIA_D_matrix(i,j) = d(4);
                RIA_D_matrix(j,i) = d(4);
            end

            idx_prof = idx_prof + 1;

            profiles{1}{idx_prof} = prof{1};    % Classic
            profiles{2}{idx_prof} = prof{2};    % RI
            profiles{3}{idx_prof} = prof{3};    % A
            profiles{4}{idx_prof} = prof{4};    % RIA
        end
    end

    %% =====================================================
    % GLOBAL DISTANCE
    %% =====================================================

    D_global = mean(D_matrix(triu(true(N),1)));
    RI_D_global = mean(RI_D_matrix(triu(true(N),1)));
    A_D_global = mean(A_D_matrix(triu(true(N),1)));
    RIA_D_global = mean(RIA_D_matrix(triu(true(N),1)));

    %% =====================================================
    % ALL TIME COORDINATES
    %% =====================================================

    %N = number_spikes*2*(factorial(num_trains)/(2*factorial(num_trains-2)));
    
    base_profile_idx = find(Distances, 1, 'first');
    t_all = zeros(sum(cellfun(@length,profiles{base_profile_idx})),1);

    idx = 1;

    if ~isempty(base_profile_idx)
        for p = 1:length(profiles{base_profile_idx})
            if ~isempty(profiles{base_profile_idx}{p})
                v = profiles{base_profile_idx}{p}(:,1);
                n = length(v);
    
                t_all(idx:idx+n-1) = v;
                idx = idx + n;
            end
        end
    end

    t_all = unique(sort(t_all));

    %% =====================================================
    % GLOBAL PROFILE SPIKE-distance
    %% =====================================================

    prof = cell(1,4);
    for i=1:4
        
        for k = 1:length(t_all)
    
            t = t_all(k);
    
            vals_left  = zeros(1,length(profiles{i}));
            vals_right = zeros(1,length(profiles{i}));
    
            has_discontinuity = false;
    
            %% -------------------------------------------------
            % scan all pairwise profiles
            %% -------------------------------------------------
    
            for p = 1:length(profiles{i})
    
                P = profiles{i}{p};

                if isempty(P)
                    continue;
                end
    
                idx = find(P(:,1) == t);
    
                %% =============================================
                % CASE 1 : discontinuity in this profile
                %% =============================================
    
                if length(idx) == 2
    
                    has_discontinuity = true;
    
                    vals_left(p)  = P(idx(1),2);
                    vals_right(p) = P(idx(2),2);
    
                %% =============================================
                % CASE 2 : single point
                %% =============================================
    
                elseif isscalar(idx) %faster than length(idx) == 1
    
                    vals_left(p)  = P(idx,2);
                    vals_right(p) = P(idx,2);
    
                %% =============================================
                % CASE 3 : interpolation
                %% =============================================
    
                else
    
                    idx_before = find(P(:,1) < t, 1, 'last');
                    idx_after  = find(P(:,1) > t, 1, 'first');
    
                    if ~isempty(idx_before) && ~isempty(idx_after)
    
                        t1 = P(idx_before,1);
                        t2 = P(idx_after,1);
    
                        S1 = P(idx_before,2);
                        S2 = P(idx_after,2);
    
                        %% avoid division by zero
                        if t2 ~= t1
                            S_interp = S1 + (S2-S1)*(t-t1)/(t2-t1);
                        else
                            S_interp = S1;
                        end
    
                        vals_left(p)  = S_interp;
                        vals_right(p) = S_interp;
    
                    end
                end
            end

            %% -------------------------------------------------
            % averaging
            %% -------------------------------------------------
    
            if has_discontinuity
                prof{i} = [ prof{i} ; t mean(vals_left) ; t mean(vals_right)];
            else
                prof{i} = [prof{i} ; t mean(vals_left)];
            end
        end
    end
    
    % sort
    profile_global = sortrows(prof{1},1);
    RI_profile_global = sortrows(prof{2},1);
    A_profile_global = sortrows(prof{3},1);
    RIA_profile_global = sortrows(prof{4},1);

end
