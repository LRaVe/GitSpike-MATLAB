
function all_shifts = f_latency_correction_sim_ann (start_spikes, tmin, tmax)

max_tau1=inf; mexy=0;

STS=SpikeTrainSet(start_spikes, tmin, tmax);
[~, ~, ~, synchProfile, ~, ~] = STS.SPIKEsynchroProfile(tmin, tmax, max_tau1);
[SPIKE_synchro_mat, ~, ~]= STS.SPIKESynchroMatrix(tmin, tmax, max_tau1);
[ssm, ~, ~, ~] = STS.Ss_So_Sto_profs(tmin, tmax, max_tau1);
num_trains=size(SPIKE_synchro_mat,1);

sim_ann_temp_fact=10000;     % 1000000 (the higher the more precision)
sim_ann_stop_diagonal=num_trains-1;          % assumes no overlap (but should be a parameter)
iter_unit=1000000;             % length of unit for piecewise dynamic allocation

start_tim =synchProfile(2,:);
%c_prof=synchProfile(1,:);
%d_prof=SorderProfile(1,:);
%e_prof=STorderProfile(1,:);
%c_value=mean(synchProfile(1,:),2);
%f_value=mean(STorderProfile(1,:),2);
ss_matches=sum(ssm,2)'/2;

all_trains = f_all_trains(start_spikes);
pairs=nchoosek(1:num_trains,2);
indies=zeros(num_trains,num_trains-1);
for trc=1:num_trains
    [indies(trc,:),~]=find(pairs==trc);
end
indies=sort(indies,2);

% in the end this should become a function with as few parameters as
% possible

% inputs / parameters:

% start_rmse, start_tim, all_trains, num_trains
% ss_matches (number of matches per spike train pair), ssm (SPIKE-synchronization matrix, number of spike train pairs x total number of spikes) 
% SPIKE_synchro_mat (Matrix of all pairwise SPIKE-synchronization values, N x N)

all_shifts=zeros(1,num_trains);

pairs=nchoosek(1:num_trains,2);
num_pairs=num_trains*(num_trains-1)/2;
            start_spike_diffs=zeros(1,num_pairs);
start_rmse=zeros(1,num_pairs);    % Random Mean Square Error
for pac=1:num_pairs
    if ss_matches(pac)>0         % to avoid start_spike_diffs(pac)==NAN
                    start_spike_diffs(pac)=mean(start_tim(ssm(pac,:) & all_trains==pairs(pac,2))-start_tim(ssm(pac,:) & all_trains==pairs(pac,1)));
        start_rmse(pac)=sqrt(mean((start_tim(ssm(pac,:) & all_trains==pairs(pac,2))-start_tim(ssm(pac,:) & all_trains==pairs(pac,1))).^2));
        % ##### here directed values, dependent on spike train order (for shifts)     % make sure that all spikes are within the bounds
    end
end
            start_spike_diffs(abs(start_spike_diffs)<1e-14)=0;
start_rmse(start_rmse<1e-14)=0;

            spike_diffs_mat=tril(ones(num_trains),-1);
            spike_diffs_mat(spike_diffs_mat==1)=start_spike_diffs;
            spike_diffs_mat=-spike_diffs_mat+spike_diffs_mat';                            % always antisymmetric !

sa_mask_mat=triu(ones(num_trains), -sim_ann_stop_diagonal) - triu(ones(num_trains), sim_ann_stop_diagonal+1);
sa_masks=sa_mask_mat(logical(tril(ones(num_trains),-1)))';
sa_start_cost = mean(start_rmse(logical(sa_masks)));

old_tim=start_tim;
old_rmse=start_rmse;
sa_old_cost=sa_start_cost;

start_rmse=start_rmse*1;     % somehow needed to preserve memory (Matlab-bug) ###
sa_costs=sa_start_cost;

keep_range=any(diag(SPIKE_synchro_mat,1)==0);                               % to check for the case where a spike train pair does not have any matches
if keep_range
    diagonal_problem=find(diag(SPIKE_synchro_mat,1)==0);
    first_indies=zeros(1,num_trains);
    last_indies=zeros(1,num_trains);
    for trc=1:num_trains
        if num_spikes(trc)>0
            first_indies(trc)=find(all_trains==trc,1,'first');
            last_indies(trc)=find(all_trains==trc,1,'last');
        end
    end
    start_separation=max(start_tim(first_indies(first_indies>0)))-min(start_tim(last_indies(last_indies>0)));  % Update in MEX: first_indies / last_indies>0
    sa_separations=start_separation;
else
    first_indies=zeros(1,num_trains);
    last_indies=zeros(1,num_trains);
end

if mexy==0 || ~exist(['SOL_sim_ann_loop_MEX.',mexext],'file')              % if MEX-version not available: Matlab-version
    sa_costs(iter_unit)=nan;
    if keep_range
        sa_separations(iter_unit)=nan;
    end
    min_cost=sa_start_cost;
    min_tim=start_tim;
    T=1;               % starting temperature (to be optimized) ######
    T_end=T/sim_ann_temp_fact;      % final temperature
    alpha=0.9;         % cooling factor
    min_iter=0;
    total_iter=1;
    sum_condi=0;
    while T>T_end
        iterations=0;
        succ_iter=0;
        
        while iterations<100*num_trains && succ_iter<10*num_trains     % ###### certain number of iterations or certain number of successful iterations
            new_tim=old_tim;
            train=randi(num_trains,1);
            displacement=randn(1)*sa_old_cost;    % random displacement (of the order of old_cost !)                                   % mean_isi/2
            
            new_tim(all_trains==train)=old_tim(all_trains==train)+displacement;
            
            new_rmse=old_rmse;
            % change only the necessary pairs: n-1 runs instead of n*(n-1)/2, factor n/2 less, which means linear order (~n) instead of squared (~n^2)
            for pac=indies(train,:)
                if ss_matches(pac)>0
                    new_rmse(pac)=sqrt(mean((new_tim(ssm(pac,:) & all_trains==pairs(pac,2))- ...
                        new_tim(ssm(pac,:) & all_trains==pairs(pac,1))).^2));
                    % ##### here directed values, dependent on spike train order (for shifts)
                end
            end
            sa_new_cost = mean(new_rmse(logical(sa_masks)));
            if keep_range
                sa_separation=max(new_tim(first_indies))-min(new_tim(last_indies));
                sa_delta_cost=sa_new_cost+10000*(sa_separation>0)-sa_old_cost;
            else
                sa_delta_cost=sa_new_cost-sa_old_cost;
            end
            
            condi=(sa_delta_cost<0 || exp(-sa_delta_cost/T)>rand(1));
            sum_condi=sum_condi+condi;
            
            if condi
                old_tim=new_tim;
                old_rmse=new_rmse;   % needed for plotting (and for correct updating!!!)
                sa_old_cost=sa_new_cost;
                succ_iter=succ_iter+1;
                if sa_new_cost<min_cost
                    min_iter=total_iter+iterations;
                    sa_min_cost=sa_new_cost;
                    min_tim=old_tim;
                end
            else
            end
            iterations=iterations+1;
            sa_costs(total_iter+iterations)=sa_old_cost;
            if keep_range
                if condi
                    sa_separations(total_iter+iterations)=sa_separation;
                else
                    sa_separations(total_iter+iterations)=sa_separations(total_iter+iterations-1);
                end
            end
        end
        total_iter=total_iter+iterations;
        if mod(total_iter,iter_unit)>iter_unit-100*num_trains
            sa_costs((ceil(total_iter/iter_unit)+1)*iter_unit)=nan;               % initialization of next 'iter_unit' iterations
            if keep_range
                sa_separations((ceil(total_iter/iter_unit)+1)*iter_unit)=nan;               % initialization of next 'iter_unit' iterations
            end
            if num_trains>10
                iteration=total_iter;
            end
        end
        T=T*alpha;   % cool down
        if succ_iter==0
            break;
        end
    end
else
    [sa_costs,total_iter,~,~,sum_condi,min_tim,sa_separations]=SOL_sim_ann_loop_MEX(int32(num_trains), int32(iter_unit), sa_old_cost, old_tim, ...
        int32(all_trains), old_rmse, int32(indies), int32(ssm), int32(ss_matches), int32(pairs), int32(sa_masks), int32(sim_ann_temp_fact), int32(keep_range), ...
        int32(first_indies), int32(last_indies));   % min_iter,min_cost
    if keep_range==1
        sa_separations(1)=start_separation;
    end
    total_iter=total_iter+1;
end
sa_costs=sa_costs(1:total_iter);
sa_last_cost=sa_costs(total_iter);
[sa_min_cost,min_iter]=min(sa_costs);

% statistics and actual shifts
min_shifts=start_tim-min_tim;
for trc=1:num_trains
    all_shifts(trc)=mean(min_shifts(all_trains==trc));
end
if any(isnan(all_shifts))                            % correction (intrapolation) for empty spike trains
    isnans=find(isnan(all_shifts));
    indies=find(~isnan(all_shifts));
    values=all_shifts(indies);
    all_shifts(isnans)=interp1(indies,values,isnans);
end
sa_end_cost=sa_min_cost;
acceptance_prob=sum_condi/total_iter*100;

