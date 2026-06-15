#include "mex.h"
#include <math.h>

/* Gateway function called directly by MATLAB */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* 1. Check input arguments */
    if (nrhs < 7) {
        mexErrMsgIdAndTxt("f_brute_force_mex:nrhs", "At least 7 inputs required (CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice).");
    }
    
    /* 2. Retrieve basic parameters from MATLAB */
    int num_neurons = (int)mxGetScalar(prhs[1]);
    int num_stimuli = (int)mxGetScalar(prhs[2]);
    int num_repetitions = (int)mxGetScalar(prhs[3]);
    double t1 = mxGetScalar(prhs[4]);
    double t2 = mxGetScalar(prhs[5]);
    
    /* Safety check against combinatorial explosion (Max 20-22 bits for integer safety) */
    if (num_neurons > 22) {
        mexErrMsgIdAndTxt("f_brute_force_mex:largeN", "N is too large for memory/integer capacity. Keep N <= 20.");
    }
    
    /* 3. Compute the total number of combinations (2^N - 1) */
    unsigned long long total_combinations = (1ULL << num_neurons) - 1;
    
    /* 4. Prepare outputs for MATLAB */
    /* Output 1: Best performance score (Scalar double) */
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *best_perf_overall = mxGetPr(plhs[1]);
    *best_perf_overall = -0.0; /* Initialized to 0 or -Inf depending on target */
    
    /* Output 2: Array to store the evaluation history */
    plhs[2] = mxCreateDoubleMatrix(1, (mwSize)total_combinations, mxREAL);
    double *history_perf_brute = mxGetPr(plhs[2]);
    
    /* 5. Pre-allocate the best mask array (internal workspace) */
    double *best_mask_overall = mxCalloc(num_neurons, sizeof(double));
    
    /* 6. Allocate an mxArray container for the current mask to pass it to the internal MATLAB calculator */
    mxArray *current_mask_mx = mxCreateDoubleMatrix(num_neurons, 1, mxREAL);
    double *current_mask_ptr = mxGetPr(current_mask_mx);
    
    /* Arrays to handle inputs/outputs when calling MATLAB's calculate_integrated_P_optimized */
    mxArray *lhs_eval[2];
    mxArray *rhs_eval[7];
    rhs_eval[0] = (mxArray*)prhs[0]; /* CellMatrix pointer */
    rhs_eval[1] = current_mask_mx;   /* Dynamic mask pointer */
    rhs_eval[2] = (mxArray*)prhs[2]; /* num_stimuli */
    rhs_eval[3] = (mxArray*)prhs[3]; /* num_repetitions */
    rhs_eval[4] = (mxArray*)prhs[4]; /* t1 */
    rhs_eval[5] = (mxArray*)prhs[5]; /* t2 */
    rhs_eval[6] = (mxArray*)prhs[6]; /* metric_choice */
    
    double max_score = -999999.0;
    
    /* 7. Main Combinatorial Loop (Binary Increment) */
    for (unsigned long long i = 1; i <= total_combinations; i++) {
        
        /* Fast bit-shifting logic: replaces dec2bin and character subtraction */
        for (int b = 0; b < num_neurons; b++) {
            /* Check if the b-th bit of integer 'i' is set to 1 */
            if ((i >> (num_neurons - 1 - b)) & 1ULL) {
                current_mask_ptr[b] = 1.0;
            } else {
                current_mask_ptr[b] = 0.0;
            }
        }
        
        /* Call MATLAB's function from inside the C code */
        mexCallMATLAB(2, lhs_eval, 7, rhs_eval, "calculate_integrated_P_optimized");
        
        /* Retrieve the computed performance score P */
        double current_perf = mxGetScalar(lhs_eval[0]);
        history_perf_brute[i - 1] = current_perf;
        
        /* Track the absolute maximum */
        if (current_perf > max_score) {
            max_score = current_perf;
            /* Save the winning pattern */
            for (int k = 0; k < num_neurons; k++) {
                best_mask_overall[k] = current_mask_ptr[k];
            }
        }
        
        /* Clean up loop iteration outputs to avoid memory leaks */
        mxDestroyArray(lhs_eval[0]);
        mxDestroyArray(lhs_eval[1]);
    }
    
    /* Set the final global max score to the output scalar */
    *best_perf_overall = max_score;
    
    /* 8. Find active neurons IDs (1-based index) from the best mask found */
    int active_count = 0;
    for (int k = 0; k < num_neurons; k++) {
        if (best_mask_overall[k] == 1.0) active_count++;
    }
    
    /* Output 0: Vector of active neuron indices (1 x active_count) */
    plhs[0] = mxCreateDoubleMatrix(1, active_count, mxREAL);
    double *best_subpop = mxGetPr(plhs[0]);
    
    int idx_out = 0;
    for (int k = 0; k < num_neurons; k++) {
        if (best_mask_overall[k] == 1.0) {
            best_subpop[idx_out] = (double)(k + 1); /* +1 because C is 0-indexed and MATLAB is 1-indexed */
            idx_out++;
        }
    }
    
    /* Free remaining persistent memory allocations */
    mxDestroyArray(current_mask_mx);
    mxFree(best_mask_overall);
}