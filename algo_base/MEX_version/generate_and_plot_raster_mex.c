#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Comparator function for quicksort (qsort) to sort spike trains chronologically */
int compare_doubles(const void *a, const void *b) {
    double temp = *(const double*)a - *(const double*)b;
    if (temp > 0) return 1;
    if (temp < 0) return -1;
    return 0;
}

/* Pseudo-random uniform generator between 0 and 1 */
double rand_double() {
    return (double)rand() / (double)RAND_MAX;
}

/* Gateway function called by MATLAB */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* 1. Input argument validation */
    if (nrhs < 8) {
        mexErrMsgIdAndTxt("generate_and_plot_raster_mex:nrhs", 
            "At least 8 arguments required (num_stimuli, num_repetitions, num_coding_neurons, num_neurons, t1, t2, base_rate, refrac).");
    }

    /* Initialize pseudo-random number generator seed */
    static int seeded = 0;
    if (!seeded) {
        srand((unsigned int)time(NULL));
        seeded = 1;
    }

    /* 2. Retrieve scalar configuration parameters */
    int num_stimuli        = (int)mxGetScalar(prhs[0]);
    int num_repetitions    = (int)mxGetScalar(prhs[1]);
    int num_coding_neurons = (int)mxGetScalar(prhs[2]);
    int num_neurons        = (int)mxGetScalar(prhs[3]);
    double t1              = mxGetScalar(prhs[4]);
    double t2              = mxGetScalar(prhs[5]);
    double base_rate       = mxGetScalar(prhs[6]);
    double refrac          = mxGetScalar(prhs[7]);

    /* 3. Create the output CellMatrix (Dimensions: num_neurons x num_stimuli x num_repetitions) */
    mwSize dims[3] = {(mwSize)num_neurons, (mwSize)num_stimuli, (mwSize)num_repetitions};
    plhs[0] = mxCreateCellArray(3, dims);

    /* Allocate temporary buffers to store generated spike timestamps */
    int max_approx_spikes = (int)((t2 - t1) * (num_coding_neurons * base_rate) * 3) + 100;
    double *pooled_buffer = (double *)mxMalloc(max_approx_spikes * sizeof(double));
    double *noise_buffer  = (double *)mxMalloc(max_approx_spikes * sizeof(double));

    /* ========================================================================= */
    /* 4. GENERATION OF CODING NEURONS (1 to c)                                  */
    /* ========================================================================= */
    double pooled_rate = (double)num_coding_neurons * base_rate;

    for (int st = 0; st < num_stimuli; st++) {
        
        /* Generate a homogeneous Poisson process with an absolute refractory period */
        double current_time = t1;
        int num_spikes = 0;

        while (current_time <= t2 && num_spikes < max_approx_spikes) {
            double u = rand_double();
            /* Avoid log(0) undefined behavior */
            if (u >= 1.0) u = 0.999999; 
            
            double interval = refrac - log(1.0 - u) / pooled_rate;
            current_time += interval;
            
            if (current_time <= t2) {
                pooled_buffer[num_spikes] = current_time;
                num_spikes++;
            }
        }

        /* Distribute the pooled spikes across each trial repetition */
        for (int rp = 0; rp < num_repetitions; rp++) {
            if (num_spikes > 0) {
                
                /* Create and shuffle indices (equivalent to MATLAB's randperm) */
                int *shuffled_indices = (int *)mxMalloc(num_spikes * sizeof(int));
                for (int i = 0; i < num_spikes; i++) shuffled_indices[i] = i;
                for (int i = num_spikes - 1; i > 0; i--) {
                    int j = rand() % (i + 1);
                    int tmp = shuffled_indices[i];
                    shuffled_indices[i] = shuffled_indices[j];
                    shuffled_indices[j] = tmp;
                }

                /* Mutually exclusive distribution among coding neurons */
                for (int nc = 0; nc < num_coding_neurons; nc++) {
                    
                    /* Count how many spikes are assigned to this specific neuron */
                    int assigned_count = 0;
                    for (int i = nc; i < num_spikes; i += num_coding_neurons) {
                        assigned_count++;
                    }

                    if (assigned_count > 0) {
                        mxArray *spike_train = mxCreateDoubleMatrix(1, assigned_count, mxREAL);
                        double *out_ptr = mxGetPr(spike_train);

                        int out_idx = 0;
                        for (int i = nc; i < num_spikes; i += num_coding_neurons) {
                            out_ptr[out_idx] = pooled_buffer[shuffled_indices[i]];
                            out_idx++;
                        }

                        /* Sort the assigned spikes chronologically */
                        qsort(out_ptr, assigned_count, sizeof(double), compare_doubles);

                        /* Calculate the linear index for MATLAB's Column-Major CellMatrix layout */
                        mwIndex cell_idx = nc + (st * num_neurons) + (rp * num_neurons * num_stimuli);
                        mxSetCell(plhs[0], cell_idx, spike_train);
                    } else {
                        mwIndex cell_idx = nc + (st * num_neurons) + (rp * num_neurons * num_stimuli);
                        mxSetCell(plhs[0], cell_idx, mxCreateDoubleMatrix(1, 0, mxREAL));
                    }
                }
                mxFree(shuffled_indices);
            } else {
                /* If no spikes were generated in the pool, populate cells with empty matrices */
                for (int nc = 0; nc < num_coding_neurons; nc++) {
                    mwIndex cell_idx = nc + (st * num_neurons) + (rp * num_neurons * num_stimuli);
                    mxSetCell(plhs[0], cell_idx, mxCreateDoubleMatrix(1, 0, mxREAL));
                }
            }
        }
    }

    /* ========================================================================= */
    /* 5. GENERATION OF NON-CODING NEURONS (c+1 to N)                            */
    /* ========================================================================= */
    for (int st = 0; st < num_stimuli; st++) {
        for (int rp = 0; rp < num_repetitions; rp++) {
            for (int nc = num_coding_neurons; nc < num_neurons; nc++) {
                
                double current_time = t1;
                int noise_spikes = 0;

                while (current_time <= t2 && noise_spikes < max_approx_spikes) {
                    double u = rand_double();
                    if (u >= 1.0) u = 0.999999;

                    double interval = refrac - log(1.0 - u) / base_rate;
                    current_time += interval;

                    if (current_time <= t2) {
                        noise_buffer[noise_spikes] = current_time;
                        noise_spikes++;
                    }
                }

                /* Export the independent noise spike train to the corresponding cell */
                mxArray *spike_train = mxCreateDoubleMatrix(1, noise_spikes, mxREAL);
                if (noise_spikes > 0) {
                    double *out_ptr = mxGetPr(spike_train);
                    memcpy(out_ptr, noise_buffer, noise_spikes * sizeof(double));
                }

                mwIndex cell_idx = nc + (st * num_neurons) + (rp * num_neurons * num_stimuli);
                mxSetCell(plhs[0], cell_idx, spike_train);
            }
        }
    }

    /* Free temporary working buffers from memory */
    mxFree(pooled_buffer);
    mxFree(noise_buffer);
}