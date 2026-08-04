#include "mex.h"
#include <math.h>

/* Gateway function called directly by MATLAB */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* 1. Check input arguments */
    if (nrhs < 7) {
        mexErrMsgIdAndTxt("f_brute_force_mex:nrhs", "At least 7 inputs required.");
    }
    
    /* 2. Retrieve basic parameters from MATLAB */
    int num_neurons = (int)mxGetScalar(prhs[1]);
    
    if (num_neurons > 20) {
        mexErrMsgIdAndTxt("f_brute_force_mex:largeN", "N is too large. Keep N <= 20.");
    }
    
    /* 3. Compute the total number of combinations (2^N - 1) */
    unsigned long long total_combinations = (1ULL << num_neurons) - 1;
    
    /* 4. Prepare outputs for MATLAB */
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *best_perf_overall = mxGetPr(plhs[1]);
    *best_perf_overall = -999999.0; 
    
    plhs[2] = mxCreateDoubleMatrix(1, (mwSize)total_combinations, mxREAL);
    double *history_perf_brute = mxGetPr(plhs[2]);
    
    /* 5. Pre-allocate the best mask array */
    double *best_mask_overall = mxCalloc(num_neurons, sizeof(double));
    
    /* 6. Allocate an mxArray container for the current mask */
    mxArray *current_mask_mx = mxCreateDoubleMatrix(num_neurons, 1, mxREAL);
    double *current_mask_ptr = mxGetPr(current_mask_mx);
    
    /* Arrays to handle inputs/outputs when calling the MEX P engine */
    mxArray *lhs_eval[1]; 
    mxArray *rhs_eval[7];
    
    /* ========================================================================= */
    /* CORRECTION CRITIQUE : On réutilise directement les pointeurs mxArray      */
    /* d'origine fournis par MATLAB (prhs). Cela garantit le maintien des types. */
    /* ========================================================================= */
    rhs_eval[0] = (mxArray*)prhs[0]; /* CellMatrix pointer */
    rhs_eval[1] = current_mask_mx;   /* Dynamic mask pointer */
    rhs_eval[2] = (mxArray*)prhs[2]; /* num_stimuli d'origine */
    rhs_eval[3] = (mxArray*)prhs[3]; /* num_repetitions d'origine */
    rhs_eval[4] = (mxArray*)prhs[4]; /* t1 d'origine */
    rhs_eval[5] = (mxArray*)prhs[5]; /* t2 d'origine */
    rhs_eval[6] = (mxArray*)prhs[6]; /* metric_choice d'origine */
    
    double max_score = -999999.0;
    
    /* 7. Main Combinatorial Loop */
    for (unsigned long long i = 1; i <= total_combinations; i++) {
        
        /* Fast bit-shifting logic */
        for (int b = 0; b < num_neurons; b++) {
            if ((i >> (num_neurons - 1 - b)) & 1ULL) {
                current_mask_ptr[b] = 1.0;
            } else {
                current_mask_ptr[b] = 0.0;
            }
        }
        
        /* Appel direct du binaire MEX */
        mexCallMATLAB(1, lhs_eval, 7, rhs_eval, "calculate_integrated_P_optimized_mex");
        
        double current_perf = mxGetScalar(lhs_eval[0]);
        history_perf_brute[i - 1] = current_perf;
        
        /* Track the absolute maximum */
        if (current_perf > max_score) {
            max_score = current_perf;
            for (int k = 0; k < num_neurons; k++) {
                best_mask_overall[k] = current_mask_ptr[k];
            }
        }
        
        /* Clean up memory */
        mxDestroyArray(lhs_eval[0]);
    }
    
    *best_perf_overall = max_score;
    
    /* 8. Find active neurons IDs (1-based index) */
    int active_count = 0;
    for (int k = 0; k < num_neurons; k++) {
        if (best_mask_overall[k] == 1.0) active_count++;
    }
    
    plhs[0] = mxCreateDoubleMatrix(1, active_count, mxREAL);
    double *best_subpop = mxGetPr(plhs[0]);
    
    int idx_out = 0;
    for (int k = 0; k < num_neurons; k++) {
        if (best_mask_overall[k] == 1.0) {
            best_subpop[idx_out] = (double)(k + 1); 
            idx_out++;
        }
    }
    
    /* Free remaining persistent memory allocations (Plus besoin de détruire les rhs_eval) */
    mxDestroyArray(current_mask_mx);
    mxFree(best_mask_overall);
}