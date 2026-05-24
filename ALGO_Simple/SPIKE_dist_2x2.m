%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026




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
            if Distances(1)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(2)==1 %RIA-SPIKE-distance
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
            if Distances(1)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes1(idx_1) S];
            end
            if Distances(2)==1 %RIA-SPIKE-distance
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
            if Distances(1)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(2)==1 %RIA-SPIKE-distance
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
            if Distances(1)==1 %A-SPIKE-distance
                S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*mean([ISI_dist_1 ISI_dist_2])*max(mean([ISI_dist_1 ISI_dist_2]),threshold));
                SPIKE_distance_profile{3}{end+1} = [spikes2(idx_2) S];
            end
            if Distances(2)==1 %RIA-SPIKE-distance
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
        end
    
        % remove duplicates
        [~, idx] = unique(profile,'rows','stable');
        profile = profile(idx,:);
    
        % sort
        profile = sortrows(profile,1);

        profile_mat{k}=profile;
    end

    

    %% =====================================================
    % FINAL DISTANCE
    %% =====================================================
    for k=1:4
        if Distances(k)==1
    
            t = profile_mat{k}(:,1);
            S = profile_mat{k}(:,2);
    
            SPIKE_distance_2x2(k) = trapz(t,S)/(t_max-t_min);
        end
    end
end



