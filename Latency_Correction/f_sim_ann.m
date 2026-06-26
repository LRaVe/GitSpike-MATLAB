function shifts = f_sim_ann(trains,tmin,tmax)
    
    n_trains = length(trains);
    shifts = zeros(1,n_trains);

    [old_matrix, old_value] = f_Cost_matrix(trains,tmin,tmax);    


    max_iter = 1000; % Maximum number of iterations
    temp = 1; % Initial temperature
    

    costs = zeros(1, max_iter); % Store costs for each iteration
    costs(1) = old_value;



