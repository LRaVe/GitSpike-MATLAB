#include "mex.h"

/*
 * C equivalent of compute_discrimination_performance.m
 * MATLAB Syntax: P = mex_compute_discrimination_performance(D, labels);
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // 1. Check for proper number of input and output arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MEX:compute_perf:InvalidInput", "Two input arguments are required (D, labels).");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MEX:compute_perf:InvalidOutput", "Only one output variable is allowed.");
    }

    // 2. Get pointers to the input data
    double *D = mxGetPr(prhs[0]);
    double *labels = mxGetPr(prhs[1]);

    // 3. Get and verify dimensions
    size_t T = mxGetM(prhs[0]); // Number of rows (Recordings)
    size_t ncols = mxGetN(prhs[0]);

    if (T != ncols) {
        mexErrMsgIdAndTxt("MEX:compute_perf:NotSquare", "The distance matrix D must be square.");
    }

    size_t len_labels = mxGetNumberOfElements(prhs[1]);
    if (len_labels != T) {
        mexErrMsgIdAndTxt("MEX:compute_perf:DimensionMismatch", "The length of 'labels' must match the dimensions of 'D'.");
    }

    // 4. Initialize accumulators
    double sum_intra = 0.0;
    double sum_inter = 0.0;
    size_t count_intra = 0;
    size_t count_inter = 0;

    // 5. Double loop to traverse the strict upper triangle (j > i)
    for (size_t i = 0; i < T; i++) {
        for (size_t j = i + 1; j < T; j++) {
            // MATLAB stores matrices in Column-Major order.
            // Element D(i+1, j+1) in MATLAB maps to D[i + j * T] in C.
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

    // 6. Compute means (safeguarded against division by zero)
    double mean_intra = (count_intra > 0) ? (sum_intra / (double)count_intra) : 0.0;
    double mean_inter = (count_inter > 0) ? (sum_inter / (double)count_inter) : 0.0;

    // 7. Create and assign the output matrix
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *P = mxGetPr(plhs[0]);
    
    // Discrimination performance
    *P = mean_inter - mean_intra;
}