#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Fonction utilitaire pour appeler calculate_integrated_P_optimized */
double call_calculate_P(const mxArray *cellMatrix, const double *mask, int num_neurons,
                        int num_stimuli, int num_repetitions, double t1, double t2, 
                        const mxArray *metric_choice) {
    
    mxArray *lhs[1];
    mxArray *rhs[7];
    mxArray *mask_array;
    double *mask_ptr;
    double p_val;
    
    /* 1. Préparation des arguments d'entrée */
    rhs[0] = (mxArray *)cellMatrix;
    
    mask_array = mxCreateDoubleMatrix(num_neurons, 1, mxREAL);
    mask_ptr = mxGetPr(mask_array);
    memcpy(mask_ptr, mask, num_neurons * sizeof(double));
    rhs[1] = mask_array;
    
    rhs[2] = mxCreateDoubleScalar((double)num_stimuli);
    rhs[3] = mxCreateDoubleScalar((double)num_repetitions);
    rhs[4] = mxCreateDoubleScalar(t1);
    rhs[5] = mxCreateDoubleScalar(t2);
    rhs[6] = (mxArray *)metric_choice;
    
    /* 2. Appel de la fonction MATLAB */
    mexCallMATLAB(1, lhs, 7, rhs, "calculate_integrated_P_optimized");
    
    /* 3. Récupération du résultat et nettoyage */
    p_val = mxGetScalar(lhs[0]);
    
    mxDestroyArray(lhs[0]);
    mxDestroyArray(rhs[1]);
    mxDestroyArray(rhs[2]);
    mxDestroyArray(rhs[3]);
    mxDestroyArray(rhs[4]);
    mxDestroyArray(rhs[5]);
    
    return p_val;
}

/* Générateur uniforme entre 0 et 1 */
double rand_uniform() {
    return (double)rand() / (double)RAND_MAX;
}

/* Cœur du fichier MEX */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Entrées */
    const mxArray *cellMatrix = prhs[0];
    int num_neurons = (int)mxGetScalar(prhs[1]);
    int num_stimuli = (int)mxGetScalar(prhs[2]);
    int num_repetitions = (int)mxGetScalar(prhs[3]);
    double t1 = mxGetScalar(prhs[4]);
    double t2 = mxGetScalar(prhs[5]);
    const mxArray *metric_choice = prhs[6];
    bool showing = mxGetLogicals(prhs[7])[0];
    
    /* Paramètres de l'algorithme */
    double cooling_factor = 0.9;        
    double alpha_threshold = 1e-5;       
    int iterations_per_temp = 5 * num_neurons; 
    int N0 = 50;
    
    /* Variables de boucle et d'état */
    double *mask_0, *temp_mask, *best_mask_overall, *next_mask;
    double p_0, best_perf_overall, temp_perf, next_perf;
    double sum_mask, theta, mean_delta, T_0;
    int count, n, iter, palier_idx, nb_iterations, unchanged_temp_cycles;
    
    /* Tableaux dynamiques pour l'historique */
    int current_max_paliers = 200;
    double *matrix_grid_history; /* Stocké à plat (Row-Major temporaire ou linéaire) */
    double *history_perf;
    
    int current_max_iter = 10000;
    double *hist_iter_P;
    double *hist_iter_bestP;
    double *hist_iter_size;
    double *hist_iter_temp;
    
    /* Initialisation du générateur pseudo-aléatoire */
    srand((unsigned int)time(NULL));
    
    /* Allocation des masques */
    mask_0 = (double *)mxCalloc(num_neurons, sizeof(double));
    temp_mask = (double *)mxCalloc(num_neurons, sizeof(double));
    best_mask_overall = (double *)mxCalloc(num_neurons, sizeof(double));
    next_mask = (double *)mxCalloc(num_neurons, sizeof(double));
    
    /* Étape 1. Initialisation aléatoire du masque */
    sum_mask = 0;
    for (int i = 0; i < num_neurons; i++) {
        mask_0[i] = (rand() % 2 == 0) ? 0.0 : 1.0;
        sum_mask += mask_0[i];
    }
    if (sum_mask == 0) mask_0[rand() % num_neurons] = 1.0;
    if (sum_mask == num_neurons) mask_0[rand() % num_neurons] = 0.0;
    
    p_0 = call_calculate_P(cellMatrix, mask_0, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice);
    
    best_perf_overall = p_0;
    memcpy(best_mask_overall, mask_0, num_neurons * sizeof(double));
    memcpy(temp_mask, mask_0, num_neurons * sizeof(double));
    temp_perf = p_0;
    
    /* Étape 2. Recherche de T_0 */
    double *delta_down = (double *)mxCalloc(N0, sizeof(double));
    count = 0;
    
    for (n = 0; n < N0; n++) {
        int idx = rand() % num_neurons;
        memcpy(next_mask, temp_mask, num_neurons * sizeof(double));
        next_mask[idx] = 1.0 - next_mask[idx];
        
        sum_mask = 0;
        for (int i = 0; i < num_neurons; i++) sum_mask += next_mask[i];
        if (sum_mask == 0 || sum_mask == num_neurons) continue;
        
        next_perf = call_calculate_P(cellMatrix, next_mask, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice);
        
        if (next_perf <= temp_perf) {
            delta_down[count] = fabs(next_perf - temp_perf);
            count++;
        }
        temp_perf = next_perf;
        memcpy(temp_mask, next_mask, num_neurons * sizeof(double));
    }
    
    mean_delta = 0.005;
    if (count > 0) {
        double sum_delta = 0;
        for (int i = 0; i < count; i++) sum_delta += delta_down[i];
        mean_delta = sum_delta / count;
    }
    T_0 = -mean_delta / log(0.95);
    mxFree(delta_down);
    
    if (T_0 <= 1e-7 || isnan(T_0)) {
        mexErrMsgIdAndTxt("f_simulated_annealing:InvalidTemperature", 
                          "The initial temperature T_0 is null, too small or NaN.");
    }
    if (showing) {
        mexPrintf("T_0 found: %.6f \n", T_0);
    }
    
    /* Allocations initiales pour l'historique */
    matrix_grid_history = (double *)mxCalloc(current_max_paliers * num_neurons, sizeof(double));
    history_perf = (double *)mxCalloc(current_max_paliers, sizeof(double));
    
    hist_iter_P = (double *)mxCalloc(current_max_iter, sizeof(double));
    hist_iter_bestP = (double *)mxCalloc(current_max_iter, sizeof(double));
    hist_iter_size = (double *)mxCalloc(current_max_iter, sizeof(double));
    hist_iter_temp = (double *)mxCalloc(current_max_iter, sizeof(double));
    
    /* Étape 3. Boucle du Recuit Simulé */
    theta = T_0;
    unchanged_temp_cycles = 0;
    palier_idx = 0;
    nb_iterations = 0;
    
    /* Reset à l'état initial après l'estimation */
    memcpy(temp_mask, mask_0, num_neurons * sizeof(double));
    temp_perf = p_0;
    
    while (theta > alpha_threshold) {
        
        /* Sécurité allocation paliers (Agrandissement dynamique si nécessaire) */
        if (palier_idx >= current_max_paliers) {
            current_max_paliers += 50;
            matrix_grid_history = (double *)mxRealloc(matrix_grid_history, current_max_paliers * num_neurons * sizeof(double));
            history_perf = (double *)mxRealloc(history_perf, current_max_paliers * sizeof(double));
        }
        
        if (showing) {
            mexPrintf("Temp: %.6f | Current P: %.4f\n", theta, temp_perf);
            mexEvalString("drawnow;"); /* Force l'affichage dans la console MATLAB */
        }
        
        for (iter = 0; iter < iterations_per_temp; iter++) {
            
            /* Sécurité allocation itérations (Agrandissement dynamique) */
            if (nb_iterations >= current_max_iter) {
                current_max_iter += 10000;
                hist_iter_P = (double *)mxRealloc(hist_iter_P, current_max_iter * sizeof(double));
                hist_iter_bestP = (double *)mxRealloc(hist_iter_bestP, current_max_iter * sizeof(double));
                hist_iter_size = (double *)mxRealloc(hist_iter_size, current_max_iter * sizeof(double));
                hist_iter_temp = (double *)mxRealloc(hist_iter_temp, current_max_iter * sizeof(double));
            }
            
            double active_count = 0;
            for (int i = 0; i < num_neurons; i++) active_count += temp_mask[i];
            memcpy(next_mask, temp_mask, num_neurons * sizeof(double));
            
            if (active_count == 1) {
                int *zero_indices = (int *)mxCalloc(num_neurons, sizeof(int));
                int z_count = 0;
                for (int i = 0; i < num_neurons; i++) {
                    if (temp_mask[i] == 0) { zero_indices[z_count] = i; z_count++; }
                }
                next_mask[zero_indices[rand() % z_count]] = 1.0;
                mxFree(zero_indices);
            } else if (active_count == num_neurons) {
                int *one_indices = (int *)mxCalloc(num_neurons, sizeof(int));
                int o_count = 0;
                for (int i = 0; i < num_neurons; i++) {
                    if (temp_mask[i] == 1) { one_indices[o_count] = i; o_count++; }
                }
                next_mask[one_indices[rand() % o_count]] = 0.0;
                mxFree(one_indices);
            } else {
                int idx_explore = rand() % num_neurons;
                next_mask[idx_explore] = 1.0 - temp_mask[idx_explore];
            }
            
            next_perf = call_calculate_P(cellMatrix, next_mask, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice);
            
            if (next_perf > temp_perf) {
                temp_perf = next_perf;
                memcpy(temp_mask, next_mask, num_neurons * sizeof(double));
            } else {
                double q = exp(-fabs(next_perf - temp_perf) / theta);
                if (rand_uniform() < q) {
                    temp_perf = next_perf;
                    memcpy(temp_mask, next_mask, num_neurons * sizeof(double));
                }
            }
            
            if (temp_perf > best_perf_overall) {
                best_perf_overall = temp_perf;
                memcpy(best_mask_overall, temp_mask, num_neurons * sizeof(double));
            }
            
            /* Enregistrement de l'historique par itération */
            hist_iter_P[nb_iterations] = temp_perf;
            hist_iter_bestP[nb_iterations] = best_perf_overall;
            
            double current_size = 0;
            for (int i = 0; i < num_neurons; i++) current_size += temp_mask[i];
            hist_iter_size[nb_iterations] = current_size;
            hist_iter_temp[nb_iterations] = theta;
            
            nb_iterations++;
        }
        
        /* Enregistrement de l'historique par palier */
        for (int i = 0; i < num_neurons; i++) {
            matrix_grid_history[palier_idx * num_neurons + i] = temp_mask[i];
        }
        history_perf[palier_idx] = temp_perf;
        palier_idx++;
        
        /* Condition d'arrêt si stagnation */
        if (palier_idx >= 2 && fabs(history_perf[palier_idx - 1] - history_perf[palier_idx - 2]) < 1e-6) {
            unchanged_temp_cycles++;
            if (unchanged_temp_cycles >= 2) {
                if (showing) {
                    mexPrintf("Exit: Performance remained unchanged for 2 consecutive temperature cycles.\n");
                }
                break;
            }
        } else {
            unchanged_temp_cycles = 0;
        }
        
        theta *= cooling_factor;
    }
    
    /* 4. Assignation des sorties vers MATLAB (Copie finale) */
    
    /* Out 0: best_mask_overall */
    plhs[0] = mxCreateDoubleMatrix(1, num_neurons, mxREAL);
    memcpy(mxGetPr(plhs[0]), best_mask_overall, num_neurons * sizeof(double));
    
    /* Out 1: best_perf_overall */
    plhs[1] = mxCreateDoubleScalar(best_perf_overall);
    
    /* Out 2: nb_iterations */
    plhs[2] = mxCreateDoubleScalar((double)nb_iterations);
    
    /* Out 3: Matrix_Grid (Attention au stockage Column-Major de MATLAB !) */
    plhs[3] = mxCreateDoubleMatrix(palier_idx, num_neurons, mxREAL);
    double *out_matrix_grid = mxGetPr(plhs[3]);
    for (int i = 0; i < palier_idx; i++) {
        for (int j = 0; j < num_neurons; j++) {
            out_matrix_grid[j * palier_idx + i] = matrix_grid_history[i * num_neurons + j];
        }
    }
    
    /* Out 4: history_perf */
    plhs[4] = mxCreateDoubleMatrix(1, palier_idx, mxREAL);
    memcpy(mxGetPr(plhs[4]), history_perf, palier_idx * sizeof(double));
    
    /* Out 5 à 8: Historiques par itération */
    plhs[5] = mxCreateDoubleMatrix(1, nb_iterations, mxREAL);
    memcpy(mxGetPr(plhs[5]), hist_iter_P, nb_iterations * sizeof(double));
    
    plhs[6] = mxCreateDoubleMatrix(1, nb_iterations, mxREAL);
    memcpy(mxGetPr(plhs[6]), hist_iter_bestP, nb_iterations * sizeof(double));
    
    plhs[7] = mxCreateDoubleMatrix(1, nb_iterations, mxREAL);
    memcpy(mxGetPr(plhs[7]), hist_iter_size, nb_iterations * sizeof(double));
    
    plhs[8] = mxCreateDoubleMatrix(1, nb_iterations, mxREAL);
    memcpy(mxGetPr(plhs[8]), hist_iter_temp, nb_iterations * sizeof(double));
    
    /* Libération de la mémoire locale */
    mxFree(mask_0);
    mxFree(temp_mask);
    mxFree(best_mask_overall);
    mxFree(next_mask);
    mxFree(matrix_grid_history);
    mxFree(history_perf);
    mxFree(hist_iter_P);
    mxFree(hist_iter_bestP);
    mxFree(hist_iter_size);
    mxFree(hist_iter_temp);
}