// #include "mex.h"
// 
// /*
//  * C equivalent of compute_discrimination_performance.m
//  * MATLAB Syntax: P = compute_discrimination_performance_mex(D, labels);
//  */
// 
// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//     // 1. Check for proper number of input and output arguments
//     if (nrhs != 2) {
//         mexErrMsgIdAndTxt("MEX:compute_perf:InvalidInput", "Two input arguments are required (D, labels).");
//     }
//     if (nlhs > 1) {
//         mexErrMsgIdAndTxt("MEX:compute_perf:InvalidOutput", "Only one output variable is allowed.");
//     }
// 
//     // 2. Get pointers to the input data
//     double *D = mxGetPr(prhs[0]);
//     double *labels = mxGetPr(prhs[1]);
// 
//     // 3. Get and verify dimensions
//     size_t T = mxGetM(prhs[0]); // Number of rows (Recordings)
//     size_t ncols = mxGetN(prhs[0]);
// 
//     if (T != ncols) {
//         mexErrMsgIdAndTxt("MEX:compute_perf:NotSquare", "The distance matrix D must be square.");
//     }
// 
//     size_t len_labels = mxGetNumberOfElements(prhs[1]);
//     if (len_labels != T) {
//         mexErrMsgIdAndTxt("MEX:compute_perf:DimensionMismatch", "The length of 'labels' must match the dimensions of 'D'.");
//     }
// 
//     // 4. Initialize accumulators
//     double sum_intra = 0.0;
//     double sum_inter = 0.0;
//     size_t count_intra = 0;
//     size_t count_inter = 0;
// 
//     // 5. Double loop to traverse the strict upper triangle (j > i)
//     for (size_t i = 0; i < T; i++) {
//         for (size_t j = i + 1; j < T; j++) {
//             // MATLAB stores matrices in Column-Major order.
//             // Element D(i+1, j+1) in MATLAB maps to D[i + j * T] in C.
//             double distance = D[i + j * T];
// 
//             if (labels[i] == labels[j]) {
//                 sum_intra += distance;
//                 count_intra++;
//             } else {
//                 sum_inter += distance;
//                 count_inter++;
//             }
//         }
//     }
// 
//     // 6. Compute means (safeguarded against division by zero)
//     double mean_intra = (count_intra > 0) ? (sum_intra / (double)count_intra) : 0.0;
//     double mean_inter = (count_inter > 0) ? (sum_inter / (double)count_inter) : 0.0;
// 
//     // 7. Create and assign the output matrix
//     plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
//     double *P = mxGetPr(plhs[0]);
// 
//     // Discrimination performance
//     *P = mean_inter - mean_intra;
// }

#include "mex.h"
#include <string.h>

/*
 * Multi-signature C-MEX implementation of compute_discrimination_performance_mex.c
 * Hardcoded default metric vector to protect unmodified MATLAB functions.
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // =========================================================================
    // SIGNATURE A: 2 Input Arguments (Direct Evaluation: D, labels)
    // =========================================================================
    if (nrhs == 2) {
        double *D = mxGetPr(prhs[0]);
        double *labels = mxGetPr(prhs[1]);
        size_t T = mxGetM(prhs[0]);

        double sum_intra = 0.0;
        double sum_inter = 0.0;
        size_t count_intra = 0;
        size_t count_inter = 0;

        for (size_t i = 0; i < T; i++) {
            for (size_t j = i + 1; j < T; j++) {
                double distance = D[i + j * T]; // Column-major indexing
                if (labels[i] == labels[j]) {
                    sum_intra += distance;
                    count_intra++;
                } else {
                    sum_inter += distance;
                    count_inter++;
                }
            }
        }

        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        double *P = mxGetPr(plhs[0]);
        *P = (count_inter > 0 ? sum_inter / (double)count_inter : 0.0) - 
             (count_intra > 0 ? sum_intra / (double)count_intra : 0.0);
        return;
    }

    // =========================================================================
    // SIGNATURE B: 7 Input Arguments (Called during search exploration loops)
    // =========================================================================
    if (nrhs < 7) {
        mexErrMsgIdAndTxt("MEX:perf:InvalidInput", "Expected 7 inputs from search loop.");
    }

    const mxArray *cellMatrix = prhs[0];
    const mxArray *mask_array = prhs[1];
    double t2 = mxGetScalar(prhs[5]);

    // --- LE STRATAGÈME : On crée un faux vecteur logique [1, 0, 0, 0] attendu par d(find(Distances)) ---
    mxArray *fake_metric = mxCreateLogicalMatrix(1, 4);
    mxLogical *fake_ptr = mxGetLogicals(fake_metric);
    fake_ptr[0] = true;   // Le find() trouvera l'indice 1
    fake_ptr[1] = false;
    fake_ptr[2] = false;
    fake_ptr[3] = false;
    // -------------------------------------------------------------------------------------------------

    double *mask_ptr = mxGetPr(mask_array);
    size_t num_neurons = mxGetNumberOfElements(mask_array);

    int active_in_mask = 0;
    for (size_t b = 0; b < num_neurons; b++) {
        if (mask_ptr[b] == 1.0) active_in_mask++;
    }

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *P = mxGetPr(plhs[0]);

    if (active_in_mask == 0) {
        *P = 0.0;
        mxDestroyArray(fake_metric);
        return;
    }

    // a. Extract active neuron indices
    mxArray *neuron_ids_array = mxCreateDoubleMatrix(1, active_in_mask, mxREAL);
    double *neuron_ids_ptr = mxGetPr(neuron_ids_array);
    int idx_c = 0;
    for (size_t b = 0; b < num_neurons; b++) {
        if (mask_ptr[b] == 1.0) {
            neuron_ids_ptr[idx_c++] = (double)(b + 1);
        }
    }

    // b. Execute: [trials, labels] = build_trials(CellMatrix, neuron_ids)
    mxArray *rhs_trials[2] = { (mxArray *)cellMatrix, neuron_ids_array };
    mxArray *lhs_trials[2];
    mexCallMATLAB(2, lhs_trials, 2, rhs_trials, "build_trials");

    // c. Execute: D = compute_population_distance_matrix(trials, t2, fake_metric, 'auto')
    mxArray *auto_str = mxCreateString("auto");
    mxArray *rhs_dist[4] = { lhs_trials[0], mxCreateDoubleScalar(t2), fake_metric, auto_str };
    mxArray *lhs_dist[1];
    mexCallMATLAB(1, lhs_dist, 4, rhs_dist, "compute_population_distance_matrix");

    // d. Fast C computation of P over the generated population distance matrix (D)
    double *D = mxGetPr(lhs_dist[0]);
    double *labels = mxGetPr(lhs_trials[1]);
    size_t T = mxGetM(lhs_dist[0]);

    double sum_intra = 0.0;
    double sum_inter = 0.0;
    size_t count_intra = 0;
    size_t count_inter = 0;

    for (size_t i = 0; i < T; i++) {
        for (size_t j = i + 1; j < T; j++) {
            double distance = D[i + j * T];
            if (labels[i] == labels[j]) {
                sum_intra += distance;
                count_intra++;
            } else {
                sum_inter += distance;
                count_inter++;
            }
        }
    }

    *P = (count_inter > 0 ? sum_inter / (double)count_inter : 0.0) - 
         (count_intra > 0 ? sum_intra / (double)count_intra : 0.0);

    /* Garbage collection */
    mxDestroyArray(neuron_ids_array);
    mxDestroyArray(lhs_trials[0]);
    mxDestroyArray(lhs_trials[1]);
    mxDestroyArray(auto_str);
    mxDestroyArray(rhs_dist[1]);
    mxDestroyArray(lhs_dist[0]);
    mxDestroyArray(fake_metric);
}