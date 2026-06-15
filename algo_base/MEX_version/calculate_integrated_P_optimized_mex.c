#include "mex.h"
#include <math.h>
#include <string.h>

/* Fonction utilitaire pour trier un tableau de doubles (QSort) requis pour l'alignement des spikes */
int compare_doubles(const void *a, const void *b) {
    double temp = *(const double*)a - *(const double*)b;
    if (temp > 0) return 1;
    if (temp < 0) return -1;
    return 0;
}

/* Fonction de recherche par dichotomie pour remplacer histcounts (simulation de l'indexation des bins) */
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
    return low + 1; /* Retourne un index 1-based pour correspondre à la logique de ton code MATLAB */
}

/* Fonction Passerelle Gateway appelée directement par MATLAB */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* 1. Vérification des arguments d'entrée */
    if (nrhs < 7) {
        mexErrMsgIdAndTxt("calculate_integrated_P_mex:nrhs", 
            "7 arguments requis (CellMatrix, selection, S, R, tmin, tmax, metric).");
    }
    
    /* 2. Récupération des paramètres de base */
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
    
    /* Allocation de la matrice de distance de sortie MatrixD (num_trials x num_trials) */
    plhs[1] = mxCreateDoubleMatrix(num_trials, num_trials, mxREAL);
    double *MatrixD = mxGetPr(plhs[1]);
    
    /* 3. Trouver les neurones sélectionnés (équivalent de find(selection == 1)) */
    int *idx_selected = mxCalloc(num_neurons, sizeof(int));
    int num_selected = 0;
    for (int n = 0; n < num_neurons; n++) {
        if (selection[n] == 1.0) {
            idx_selected[num_selected] = n;
            num_selected++;
        }
    }
    
    /* Si aucun neurone n'est sélectionné, P = -Inf */
    if (num_selected == 0) {
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(plhs[0]) = -INFINITY;
        mxFree(idx_selected);
        return;
    }
    
    /* 4. Pré-extraction et fusion des trains de spikes (Summed Population) */
    /* On crée des listes dynamiques pour stocker le train fusionné de chaque essai */
    double **Precomputed_Trains = mxCalloc(num_trials, sizeof(double*));
    int *Train_Lengths = mxCalloc(num_trials, sizeof(int));
    
    for (int t = 0; t < num_trials; t++) {
        int st = t / R;
        int rp = t % R;
        
        /* Calcul de la taille totale requise pour cet essai */
        int total_spikes_est = 0;
        for (int n = 0; n < num_selected; n++) {
            int neuron_idx = idx_selected[n];
            /* Calcul de l'index linéaire dans la CellMatrix à 3 dimensions de MATLAB */
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
            /* Tri obligatoire du train fusionné (équivalent de sort() en MATLAB) */
            qsort(Precomputed_Trains[t], total_spikes_est, sizeof(double), compare_doubles);
        }
    }
    
    /* 5. Double boucle imbriquée pour le calcul des distances */
    for (int t_a = 0; t_a < num_trials; t_a++) {
        double *train_A = Precomputed_Trains[t_a];
        int len_A = Train_Lengths[t_a];
        
        for (int t_b = t_a + 1; t_b < num_trials; t_b++) {
            double *train_B = Precomputed_Trains[t_b];
            int len_B = Train_Lengths[t_b];
            
            double dval = 0.5; // Sécurité par défaut
            
            if (len_A == 0 && len_B == 0) {
                dval = 0.0;
            } else if (len_A == 0 || len_B == 0) {
                dval = 1.0;
            } else {
                if (strcasecmp(metric, "SPIKE_DISTANCE") == 0) {
                    /* --- Implémentation C native ultra-rapide de SPIKE_DISTANCE --- */
                    
                    /* Construction de la grille unifiée des temps uniques */
                    int max_edges = len_A + len_B + 2;
                    double *t_all = mxCalloc(max_edges, sizeof(double));
                    t_all[0] = tmin;
                    int edge_count = 1;
                    
                    /* Fusion et suppression des doublons simplifiée */
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
                    
                    /* Définition des intervalles (t_diff) et milieux (t_mids) */
                    int len_mids = edge_count - 1;
                    double total_integrated = 0.0;
                    
                    /* Construction des bords complets pour histcounts */
                    double *edges_A = mxCalloc(len_A + 2, sizeof(double));
                    edges_A[0] = -INFINITY; memcpy(&edges_A[1], train_A, len_A * sizeof(double)); edges_A[len_A + 1] = INFINITY;
                    
                    double *edges_B = mxCalloc(len_B + 2, sizeof(double));
                    edges_B[0] = -INFINITY; memcpy(&edges_B[1], train_B, len_B * sizeof(double)); edges_B[len_B + 1] = INFINITY;
                    
                    for (int k = 0; k < len_mids; k++) {
                        double t_diff = t_all[k+1] - t_all[k];
                        double t_mid = (t_all[k] + t_all[k+1]) / 2.0;
                        
                        /* Binning pour train A */
                        int bin_A = find_bin_index(t_mid, edges_A, len_A + 2);
                        int idx_p_A = bin_A - 1; if (idx_p_A < 1) idx_p_A = 1;
                        int idx_n_A = bin_A;     if (idx_n_A > len_A) idx_n_A = len_A;
                        
                        double x_p = train_A[idx_p_A - 1]; if (bin_A - 1 < 1) x_p = tmin;
                        double x_a = train_A[idx_n_A - 1]; if (bin_A > len_A) x_a = tmax;
                        double isi_x = x_a - x_p;
                        
                        /* Binning pour train B */
                        int bin_B = find_bin_index(t_mid, edges_B, len_B + 2);
                        int idx_p_B = bin_B - 1; if (idx_p_B < 1) idx_p_B = 1;
                        int idx_n_B = bin_B;     if (idx_n_B > len_B) idx_n_B = len_B;
                        
                        double y_p = train_B[idx_p_B - 1]; if (bin_B - 1 < 1) y_p = tmin;
                        double y_a = train_B[idx_n_B - 1]; if (bin_B > len_B) y_a = tmax;
                        double isi_y = y_a - y_p;
                        
                        /* Recherche des cibles (Spikes les plus proches) */
                        double target_x = ((t_mid - x_p) < (x_a - t_mid)) ? x_p : x_a;
                        double target_y = ((t_mid - y_p) < (y_a - t_mid)) ? y_p : y_a;
                        
                        /* Plus proche voisin de target_x dans train_B */
                        int b_bin = find_bin_index(target_x, edges_B, len_B + 2);
                        int b_p = b_bin - 1; if (b_p < 1) b_p = 1;
                        int b_n = b_bin;     if (b_n > len_B) b_n = len_B;
                        double near_B = train_B[b_n - 1];
                        if ((b_bin - 1 >= 1) && (fabs(target_x - train_B[b_p - 1]) < fabs(train_B[b_n - 1] - target_x))) {
                            near_B = train_B[b_p - 1];
                        }
                        double min_dxy = fabs(target_x - near_B);
                        
                        /* Plus proche voisin de target_y dans train_A */
                        int a_bin = find_bin_index(target_y, edges_A, len_A + 2);
                        int a_p = a_bin - 1; if (a_p < 1) a_p = 1;
                        int a_n = a_bin;     if (a_n > len_A) a_n = len_A;
                        double near_A = train_A[a_n - 1];
                        if ((a_bin - 1 >= 1) && (fabs(target_y - train_A[a_p - 1]) < fabs(train_A[a_n - 1] - target_y))) {
                            near_A = train_A[a_p - 1];
                        }
                        double min_dyx = fabs(target_y - near_A);
                        
                        /* Profils de distance */
                        double S_x = ((t_mid - x_p) * min_dyx + (x_a - t_mid) * min_dyx) / isi_x;
                        double S_y = ((t_mid - y_p) * min_dxy + (y_a - t_mid) * min_dxy) / isi_y;
                        double S_t = (S_x * isi_y + S_y * isi_x) / ((isi_x + isi_y) * (isi_x > isi_y ? isi_x : isi_y));
                        
                        total_integrated += S_t * t_diff;
                    }
                    
                    dval = total_integrated / (tmax - tmin);
                    
                    mxFree(t_all);
                    mxFree(edges_A);
                    mxFree(edges_B);
                }
            }
            
            /* Remplissage symétrique de la matrice MatrixD */
            MatrixD[t_a + t_b * num_trials] = dval;
            MatrixD[t_b + t_a * num_trials] = dval;
        }
    }
    
    /* 6. Calcul final du score de performance globale P */
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
    
    /* Sortie 0 : Score P final scalar */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = mean_inter - mean_intra;
    
    /* Nettoyage de la mémoire persistante locale */
    for (int t = 0; t < num_trials; t++) {
        if (Train_Lengths[t] > 0) mxFree(Precomputed_Trains[t]);
    }
    mxFree(Precomputed_Trains);
    mxFree(Train_Lengths);
    mxFree(idx_selected);
}