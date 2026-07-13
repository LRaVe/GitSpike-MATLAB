#include "mex.h"
#include <math.h>
#include <string.h>

#ifdef _WIN32
#define strcasecmp _stricmp
#endif

/* Utility function to sort an array of doubles (QSort) required for spike alignment */
int compare_doubles(const void *a, const void *b) {
    double temp = *(const double*)a - *(const double*)b;
    if (temp > 0) return 1;
    if (temp < 0) return -1;
    return 0;
}

/* Binary search function replacing histcounts (simulates bin indexing) */
int find_bin_index(double value, const double *edges, int num_edges) {
    if (value < edges[0]) return 1;
    if (value >= edges[num_edges - 1]) return num_edges;
    
    int low = 0, high = num_edges - 1;
    while (high - low > 1) {
        int mid = (low + high) / 2;
        if (value >= edges[mid]) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return low + 1; /* Returns a 1-based index to match MATLAB's logical indexing */
}

/* Helper function to compute the instantaneous SPIKE-distance profile at a single timestamp 't' */
double compute_instantaneous_S(double t, const double *train_A, int len_A, const double *edges_A,
                               const double *train_B, int len_B, const double *edges_B, double tmin, double tmax) {
    
    % --- Train A Profiling ---
    int bin_A = find_bin_index(t, edges_A, len_A + 2);
    int idx_p_A = bin_A - 1; if (idx_p_A < 1) idx_p_A = 1;
    int idx_n_A = bin_A;     if (idx_n_A > len_A) idx_n_A = len_A;
    
    double x_p = train_A[idx_p_A - 1]; if (bin_A - 1 < 1) x_p = tmin;
    double x_a = train_A[idx_n_A - 1]; if (bin_A > len_A) x_a = tmax;
    double isi_x = x_a - x_p;
    
    % Distance to the nearest spike in Train A
    double min_val_A = fabs(t - train_A[0]);
    for (int i = 1; i < len_A; i++) {
        double d = fabs(t - train_A[i]);
        if (d < min_val_A) min_val_A = d;
    }
    
    % --- Train B Profiling ---
    int bin_B = find_bin_index(t, edges_B, len_B + 2);
    int idx_p_B = bin_B - 1; if (idx_p_B < 1) idx_p_B = 1;
    int idx_n_B = bin_B;     if (idx_n_B > len_B) idx_n_B = len_B;
    
    double y_p = train_B[idx_p_B - 1]; if (bin_B - 1 < 1) y_p = tmin;
    double y_a = train_B[idx_n_B - 1]; if (bin_B > len_B) y_a = tmax;
    double isi_y = y_a - y_p;
    
    % Distance to the nearest spike in Train B
    double min_val_B = fabs(t - train_B[0]);
    for (int i = 1; i < len_B; i++) {
        double d = fabs(t - train_B[i]);
        if (d < min_val_B) min_val_B = d;
    }
    
    % --- Local Kreuz Distances ---
    double S_x = min_val_B; % Simplified from: ((t - x_p) * min_val_B + (x_a - t) * min_val_B) / isi_x
    double S_y = min_val_A; % Simplified from: ((t - y_p) * min_val_A + (y_a - t) * min_val_A) / isi_y
    
    % Final bivariate profile value at timestamp t
    double S_t = (S_x * isi_y + S_y * isi_x) / ((isi_x + isi_y) * (isi_x > isi_y ? isi_x : isi_y));
    return S_t;
}

/* Gateway Mex function called by MATLAB */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* 1. Input argument verification */
    if (nrhs < 7) {
        mexErrMsgIdAndTxt("calculate_integrated_P_mex:nrhs", 
            "7 arguments required (CellMatrix, selection, S, R, tmin, tmax, metric).");
    }
    
    /* 2. Retrieve base input variables */
    const mxArray *CellMatrix = prhs[0];
    double *selection = mxGetPr(prhs[1]);
    int S = (int)mxGetScalar(prhs[2]);
    int R = (int)mxGetScalar(prhs[3]);
    double tmin = mxGetScalar(prhs[4]);
    double tmax = mxGetScalar(prhs[5]);
    
    char metric[64];
    mxGetString(prhs[6], metric, sizeof(metric));
    
    int num_neurons = mxGetM(CellMatrix);
    int num_trials = S * R;
    
    /* Allocate output distance matrix MatrixD (num_trials x num_trials) */
    plhs[1] = mxCreateDoubleMatrix(num_trials, num_trials, mxREAL);
    double *MatrixD = mxGetPr(plhs[1]);
    
    /* 3. Find selected neurons (Equivalent to find(selection == 1)) */
    int *idx_selected = mxCalloc(num_neurons, sizeof(int));
    int num_selected = 0;
    for (int n = 0; n < num_neurons; n++) {
        if (selection[n] == 1.0) {
            idx_selected[num_selected] = n;
            num_selected++;
        }
    }
    
    /* If no neuron is selected, return P = -Inf */
    if (num_selected == 0) {
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(plhs[0]) = -INFINITY;
        mxFree(idx_selected);
        return;
    }
    
    /* 4. Pre-extraction and merging of spike trains (Summed Population) */
    double **Precomputed_Trains = mxCalloc(num_trials, sizeof(double*));
    int *Train_Lengths = mxCalloc(num_trials, sizeof(int));
    
    for (int t = 0; t < num_trials; t++) {
        int st = t / R;
        int rp = t % R;
        
        int total_spikes_est = 0;
        for (int n = 0; n < num_selected; n++) {
            int neuron_idx = idx_selected[n];
            mwIndex cell_subs[3] = {neuron_idx, st, rp};
            mwSize cell_index = mxCalcSingleSubscript(CellMatrix, 3, cell_subs);
            mxArray *cell_element = mxGetCell(CellMatrix, cell_index);
            if (cell_element && !mxIsEmpty(cell_element)) {
                total_spikes_est += mxGetNumberOfElements(cell_element);
            }
        }
        
        Train_Lengths[t] = total_spikes_est;
        if (total_spikes_est > 0) {
            Precomputed_Trains[t] = mxCalloc(total_spikes_est, sizeof(double));
            int current_offset = 0;
            
            for (int n = 0; n < num_selected; n++) {
                int neuron_idx = idx_selected[n];
                mwIndex cell_subs[3] = {neuron_idx, st, rp};
                mwSize cell_index = mxCalcSingleSubscript(CellMatrix, 3, cell_subs);
                mxArray *cell_element = mxGetCell(CellMatrix, cell_index);
                if (cell_element && !mxIsEmpty(cell_element)) {
                    int elements = mxGetNumberOfElements(cell_element);
                    double *spikes_ptr = mxGetPr(cell_element);
                    memcpy(&Precomputed_Trains[t][current_offset], spikes_ptr, elements * sizeof(double));
                    current_offset += elements;
                }
            }
            /* Sort combined population trains */
            qsort(Precomputed_Trains[t], total_spikes_est, sizeof(double), compare_doubles);
        }
    }
    
    /* 5. Nested loops for pairwise distance calculation */
    for (int t_a = 0; t_a < num_trials; t_a++) {
        double *train_A = Precomputed_Trains[t_a];
        int len_A = Train_Lengths[t_a];
        
        for (int t_b = t_a + 1; t_b < num_trials; t_b++) {
            double *train_B = Precomputed_Trains[t_b];
            int len_B = Train_Lengths[t_b];
            
            double dval = 0.5; 
            
            if (len_A == 0 && len_B == 0) {
                dval = 0.0;
            } else if (len_A == 0 || len_B == 0) {
                dval = 1.0;
            } else {
                if (strcasecmp(metric, "SPIKE_DISTANCE") == 0) {
                    
                    /* Build unified time grid of all critical unique points */
                    int max_edges = len_A + len_B + 2;
                    double *t_all = mxCalloc(max_edges, sizeof(double));
                    t_all[0] = tmin;
                    int edge_count = 1;
                    
                    int idxA = 0, idxB = 0;
                    while (idxA < len_A || idxB < len_B) {
                        double next_val;
                        if (idxA < len_A && idxB < len_B) {
                            if (train_A[idxA] < train_B[idxB]) { next_val = train_A[idxA++]; }
                            else if (train_A[idxA] > train_B[idxB]) { next_val = train_B[idxB++]; }
                            else { next_val = train_A[idxA++]; idxB++; }
                        } else if (idxA < len_A) {
                            next_val = train_A[idxA++];
                        } else {
                            next_val = train_B[idxB++];
                        }
                        if (next_val > tmin && next_val < tmax) {
                            t_all[edge_count++] = next_val;
                        }
                    }
                    t_all[edge_count++] = tmax;
                    
                    /* Construct bounding edges for interval tracking */
                    double *edges_A = mxCalloc(len_A + 2, sizeof(double));
                    edges_A[0] = -INFINITY; memcpy(&edges_A[1], train_A, len_A * sizeof(double)); edges_A[len_A + 1] = INFINITY;
                    
                    double *edges_B = mxCalloc(len_B + 2, sizeof(double));
                    edges_B[0] = -INFINITY; memcpy(&edges_B[1], train_B, len_B * sizeof(double)); edges_B[len_B + 1] = INFINITY;
                    
                    /* Exact Continuous Integration via Trapezoidal Method */
                    double total_integrated = 0.0;
                    
                    % Pre-calculate the instantaneous profile value at the first point
                    double S_current = compute_instantaneous_S(t_all[0], train_A, len_A, edges_A, train_B, len_B, edges_B, tmin, tmax);
                    
                    for (int k = 0; k < edge_count - 1; k++) {
                        double t_start = t_all[k];
                        double t_end = t_all[k+1];
                        double t_diff = t_end - t_start;
                        
                        % Compute the profile value at the end of the current segment
                        double S_next = compute_instantaneous_S(t_end, train_A, len_A, edges_A, train_B, len_B, edges_B, tmin, tmax);
                        
                        % Trapezoidal integration addition
                        total_integrated += ((S_current + S_next) / 2.0) * t_diff;
                        
                        % Move window forward
                        S_current = S_next;
                    }
                    
                    dval = total_integrated / (tmax - tmin);
                    
                    mxFree(t_all);
                    mxFree(edges_A);
                    mxFree(edges_B);
                }
            }
            
            /* Symmetrically fill distance matrix */
            MatrixD[t_a + t_b * num_trials] = dval;
            MatrixD[t_b + t_a * num_trials] = dval;
        }
    }
    
    /* 6. Calculate global performance score P */
    double sum_inter = 0.0, sum_intra = 0.0;
    int count_inter = 0, count_intra = 0;
    
    for (int t_a = 0; t_a < num_trials; t_a++) {
        int stim_A = t_a / R;
        for (int t_b = t_a + 1; t_b < num_trials; t_b++) {
            int stim_B = t_b / R;
            double dval = MatrixD[t_a + t_b * num_trials];
            
            if (stim_A == stim_B) {
                sum_intra += dval;
                count_intra++;
            } else {
                sum_inter += dval;
                count_inter++;
            }
        }
    }
    
    double mean_inter = (count_inter > 0) ? (sum_inter / count_inter) : 0.0;
    double mean_intra = (count_intra > 0) ? (sum_intra / count_intra) : 0.0;
    
    /* Output 0: Scalar final P score */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = mean_inter - mean_intra;
    
    /* Cleanup allocation memory */
    for (int t = 0; t < num_trials; t++) {
        if (Train_Lengths[t] > 0) mxFree(Precomputed_Trains[t]);
    }
    mxFree(Precomputed_Trains);
    mxFree(Train_Lengths);
    mxFree(idx_selected);
}