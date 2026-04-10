#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mt19937ar.h"
#include "MemoryOperation.h"
#include "include/stats.hpp"
#include <cmath>
#include <random>
#include <omp.h>
#include <vector>
#include <algorithm>
// Thread-safe Laplace noise generator
double LaplaceNoiseThreadSafe(double b, std::mt19937_64& local_rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(local_rng) - 0.5;
    double sgn = (u > 0.0) ? 1.0 : -1.0;
    
    double val = 1.0 - 2.0 * std::fabs(u);
    if (val <= 0.0) val = 1e-10; 
    return -b * sgn * std::log(val);
}
// Helper function to generate Laplace noise
// genrand_real2() is assumed to generate a uniform random double in [0, 1)
double LaplaceNoise(double b) {
    double u = genrand_real2() - 0.5;
    double sgn = (u > 0.0) ? 1.0 : -1.0;
    // Prevent log(0) in edge cases
    double val = 1.0 - 2.0 * std::fabs(u);
    if (val <= 0.0) val = 1e-10; 
    return -b * sgn * std::log(val);
}
using namespace std;

#define MeasureTime	1

string EdgeFile;
int NodeNum;
double EpsT;	// epsilon
double EpsD;	// epsilon for degree (Alg = 4)
double Delta;	// delta
string EpsT_s, Delta_s;
int PairNum;
int ItrNum;
int Alg;
double AlgPrm;
int NumericalBound;
string Alg_s;
int Bip;

// Parameters in the 2-rounds local algorithm [Imola+, USENIX22]
double TClip;
double EClip;
double EpsNsDeg;
double Eps1st;
double Eps2ndTrSt;
double Mu;

// Initialization of statslib
stats::rand_engine_t engine(1776);

FILE *FileOpen(string filename, const char *mode) {
	FILE *fp;

	if ((fp = fopen(filename.c_str(), mode)) == NULL) {
		cout << "cannot open " << filename << endl;
		exit(-1);
	}
	return fp;
}

bool checkFileExistence(const std::string& str) {
    std::ifstream ifs(str);
    return ifs.is_open();
}

// Randomly generate 0, 1, 2, ..., size-1, and store the first num values into rndperm
void MakeRndPerm(int *rndperm, int size, int num) {
	int rnd;
	int *ordperm;
	int i, j;

	// 0, 1, 2, ..., size-1 --> ordperm
	ordperm = (int *)malloc(size * sizeof(int));
	for (i = 0; i < size; i++) {
		ordperm[i] = i;
	}

	for (i = 0; i < num; i++) {
		rnd = genrand_int32() % (size - i);
		rndperm[i] = ordperm[rnd];
		for (j = rnd + 1; j < size - i; j++) {
			ordperm[j - 1] = ordperm[j];
		}
	}

	free(ordperm);
}

// Read edges from the edge file and populate both map and vector representations
void ReadEdges(map<int, int> *a_mat, vector<int> *adj_list, int *node_order){
    int node1, node2;
    int i;
    char s[1025];
    char *tok;
    FILE *fp;
    int type1, type2;

    fp = FileOpen(EdgeFile, "r");
    // Skip the first 3 lines (header)
    for(i=0; i<3; i++) fgets(s, 1024, fp);
    
    while(fgets(s, 1024, fp) != NULL){
        // 1st node --> node1
        tok = strtok(s, ",");
        node1 = atoi(tok);
        // 2nd node --> node2
        tok = strtok(NULL, ",");
        node2 = atoi(tok);
        if(node1 == node2) continue;

        // Make a bipartite graph (type = 1 (odd) or 0 (even))
        if(Bip == 1){
            type1 = node_order[node1] % 2;
            type2 = node_order[node2] % 2;
            if(type1 == type2) continue;
        }

        // If both nodes exist, add the edge
        if(node_order[node1] < NodeNum && node_order[node2] < NodeNum){
            int u = node_order[node1];
            int v = node_order[node2];
            
            // 1. Map for DP algorithms (Kept original for compatibility)
            a_mat[u][v] = 1;
            a_mat[v][u] = 1;
            
            // 2. Vector for fast true graph computations
            adj_list[u].push_back(v);
            adj_list[v].push_back(u);
        }
    }
    fclose(fp);

    // Sort and remove duplicates in adj_list to perfectly match map's behavior
    for(i = 0; i < NodeNum; i++) {
        sort(adj_list[i].begin(), adj_list[i].end());
        adj_list[i].erase(unique(adj_list[i].begin(), adj_list[i].end()), adj_list[i].end());
    }
}
// Calculate #4-cycles in the shuffle model
double CalImolaCy4(map<int, int> *a_mat){
	int pair;
	int *node_order;
	int node1, node2;
	double p1w, q1w;
	int wedge, wedge_ns;
	double rnd;
	int i, j;
	double *a_mat_node1, *a_mat_node2;
	map<int, int>::iterator aitr;
	int pair_num;
	int wedge_num;
	double c1;
	double term1, term2;

	// Initialization
	double cy4_num_ns = 0.0;

	// malloc
	malloc1D(&node_order, NodeNum);
	malloc1D(&a_mat_node1, NodeNum);
	malloc1D(&a_mat_node2, NodeNum);

	// local model
	// Flip probability (wedge) --> q1w
	q1w = 1.0 / (exp(EpsT) + 1.0);
	p1w = 1 - q1w;

	// Randomly generate 0, 1, 2, ..., NodeNum-1
	MakeRndPerm(node_order, NodeNum, NodeNum);

	// For each pair of two nodes
	pair = 0;
	for(i=0;i<NodeNum;i+=2){
		// Initialization
		wedge_num = 0;

		// Two nodes --> node1, node2
		node1 = node_order[i];
		node2 = node_order[i+1];

		// a_mat[node1] --> a_mat_node1
		for(j=0;j<NodeNum;j++) a_mat_node1[j] = 0;
		for (aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
			a_mat_node1[aitr->first] = 1;
		}
		// a_mat[node2] --> a_mat_node2
		for(j=0;j<NodeNum;j++) a_mat_node2[j] = 0;
		for (aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
			a_mat_node2[aitr->first] = 1;
		}

		// For each node
		c1 = 0.0;
		for(j=0;j<NodeNum;j++){
			if(j == node1 || j == node2) continue;

			// Original wedge --> wedge
			if (a_mat_node1[j] == 1 && a_mat_node2[j] == 1) wedge = 1;
			else wedge = 0;

			// Noisy wedge --> wedge_ns
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < q1w && wedge == 0) wedge_ns = 1;
			// 1 --> 1 (not flip)
			else if(rnd >= q1w && wedge == 1) wedge_ns = 1;
			else wedge_ns = 0;

			// Update the unbiased estimate of #wedges --> c1
			c1 += ((double)wedge_ns - (1.0 - p1w)) / (2.0 * p1w - 1.0);
		}

		// Update the unbiased estimate of #4-cycles --> cy4_num_ns
		term1 = (double)c1 * ((double)c1 - 1.0) / 2.0;
		term2 = (((double)NodeNum - 2.0) * p1w * (1.0 - p1w)) / (2.0 * (2.0 * p1w - 1.0) * (2.0 * p1w - 1.0));
		cy4_num_ns += (term1 - term2);

		pair++;
		if(PairNum != -1 && PairNum == pair) break;
	}

	// Number of pairs --> pairnum
	if(PairNum != -1) pair_num = PairNum;
	else pair_num = (int)(NodeNum / 2);

	// Update the unbiased estimate of #4-cycles --> cy4_num_ns
	

	free1D(node_order);
	free1D(a_mat_node1);
	free1D(a_mat_node2);
	return ((double)NodeNum * ((double)NodeNum - 1.0) / (4.0 * (double)pair_num)) * cy4_num_ns;
}
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <map>
#include <omp.h>

extern int NodeNum;
extern double EpsT;
// helper function: based on 6-bit mask, return the respective ID (0-10)
int get_motif_from_mask(int mask) {
    int deg[4] = {0, 0, 0, 0};
    if (mask & (1 << 0)) { deg[0]++; deg[1]++; } // e0: (0,1)
    if (mask & (1 << 1)) { deg[0]++; deg[2]++; } // e1: (0,2)
    if (mask & (1 << 2)) { deg[0]++; deg[3]++; } // e2: (0,3)
    if (mask & (1 << 3)) { deg[1]++; deg[2]++; } // e3: (1,2)
    if (mask & (1 << 4)) { deg[1]++; deg[3]++; } // e4: (1,3)
    if (mask & (1 << 5)) { deg[2]++; deg[3]++; } // e5: (2,3)
    
    std::sort(deg, deg + 4);
    
    if (deg[0]==0 && deg[1]==0 && deg[2]==0 && deg[3]==0) return 0;  // 0
    if (deg[0]==0 && deg[1]==0 && deg[2]==1 && deg[3]==1) return 1;  // 1
    if (deg[0]==1 && deg[1]==1 && deg[2]==1 && deg[3]==1) return 2;  // 2
    if (deg[0]==0 && deg[1]==1 && deg[2]==1 && deg[3]==2) return 3;  // 2 (path)
    if (deg[0]==1 && deg[1]==1 && deg[2]==1 && deg[3]==3) return 4;  // 3 (star)
    if (deg[0]==1 && deg[1]==1 && deg[2]==2 && deg[3]==2) return 5;  // 3 (3-path)
    if (deg[0]==0 && deg[1]==2 && deg[2]==2 && deg[3]==2) return 6;  // 3 (triangle)
    if (deg[0]==1 && deg[1]==2 && deg[2]==2 && deg[3]==3) return 7;  // 4 (Paw)
    if (deg[0]==2 && deg[1]==2 && deg[2]==3 && deg[3]==3) return 8;  // 5 (Diamond)
    if (deg[0]==3 && deg[1]==3 && deg[2]==3 && deg[3]==3) return 9;  // 6 (Clique)
    if (deg[0]==2 && deg[1]==2 && deg[2]==2 && deg[3]==2) return 10; // 4 (4-Cycle)
    
    return -1;
}

// helper function: count number of bits
int count_bits(int n) {
    int count = 0;
    while (n) {
        count += n & 1;
        n >>= 1;
    }
    return count;
}

// Gaussian elimination to solve A * x = b
std::vector<double> solve_linear_system(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                max_row = k;
            }
        }
        std::swap(A[i], A[max_row]);
        std::swap(b[i], b[max_row]);
        double pivot = A[i][i];
        if (std::abs(pivot) < 1e-12) continue;
        for (int j = i; j < n; j++) A[i][j] /= pivot;
        b[i] /= pivot;
        
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = i; j < n; j++) A[k][j] -= factor * A[i][j];
                b[k] -= factor * b[i];
            }
        }
    }
    return b;
}

// ---------------------------------------------------------
// One-Round Motif-Transformation Baseline (OneR)
// ---------------------------------------------------------
double CalOneR(std::map<int, int> *a_mat) {
    // Randomized Response
    double p = 1.0 / (exp(EpsT) + 1.0); 

    // adjacency matrix for generating G'
    std::vector<std::vector<char>> G_prime(NodeNum, std::vector<char>(NodeNum, 0));

    #pragma omp parallel
    {
        std::mt19937_64 local_rng(1776 + omp_get_thread_num());
        std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < NodeNum; i++) {
            for (int j = i + 1; j < NodeNum; j++) {
                int true_edge = (a_mat[i].find(j) != a_mat[i].end()) ? 1 : 0;
                double rnd = uniform_dist(local_rng);
                
                int noisy_edge = 0;
                if (true_edge == 1 && rnd >= p) noisy_edge = 1;
                else if (true_edge == 0 && rnd < p) noisy_edge = 1;
                
                G_prime[i][j] = noisy_edge;
                G_prime[j][i] = noisy_edge;
            }
        }
    }

    // enumerate subgraphs in G'
    std::vector<double> n_prime(11, 0.0);
    
    #pragma omp parallel
    {
        std::vector<double> local_n_prime(11, 0.0);

        #pragma omp for schedule(dynamic, 10)
        for (int i = 0; i < NodeNum; i++) {
            for (int j = i + 1; j < NodeNum; j++) {
                for (int k = j + 1; k < NodeNum; k++) {
                    for (int l = k + 1; l < NodeNum; l++) {
                        // sample 6 edges to form bitmask
                        int mask = 0;
                        if (G_prime[i][j]) mask |= (1 << 0);
                        if (G_prime[i][k]) mask |= (1 << 1);
                        if (G_prime[i][l]) mask |= (1 << 2);
                        if (G_prime[j][k]) mask |= (1 << 3);
                        if (G_prime[j][l]) mask |= (1 << 4);
                        if (G_prime[k][l]) mask |= (1 << 5);
                        
                        int m = get_motif_from_mask(mask);
                        local_n_prime[m] += 1.0;
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            for (int m = 0; m < 11; m++) {
                n_prime[m] += local_n_prime[m];
            }
        }
    }

    // Construct transition matrix
    std::vector<std::vector<double>> T(11, std::vector<double>(11, 0.0));
    
    // canonical mask
    int canonical_mask[11];
    std::fill_n(canonical_mask, 11, -1);
    for (int mask = 0; mask < 64; ++mask) {
        int m = get_motif_from_mask(mask);
        if (canonical_mask[m] == -1) canonical_mask[m] = mask;
    }

    for (int i = 0; i < 11; i++) {
        int true_mask = canonical_mask[i];
        for (int obs_mask = 0; obs_mask < 64; ++obs_mask) {
            int obs_m = get_motif_from_mask(obs_mask);
            int flips = count_bits(true_mask ^ obs_mask); // edges needed to flip
            double prob = std::pow(p, flips) * std::pow(1.0 - p, 6 - flips);
            T[i][obs_m] += prob;
        }
    }

    // Calculate \hat{n} = n' * T^{-1}
    std::vector<std::vector<double>> TT(11, std::vector<double>(11, 0.0));
    for (int i = 0; i < 11; i++) {
        for (int j = 0; j < 11; j++) {
            TT[i][j] = T[j][i];
        }
    }
    std::vector<double> n_hat = solve_linear_system(TT, n_prime);

    // return unbiased estimation
    double final_estimate = n_hat[10];

    if (final_estimate < 0.0) {
        final_estimate = 0.0;
    }
    
    return final_estimate;
}
double CalCN4C(map<int, int> *a_mat) {
    int *node_order;
    
    // Privacy parameters
    double eps0, eps1, eps2;
    double p, V0, gamma, lap_var_Y;
    
    // Initialization
    double cy4_num_ns = 0.0;

    // malloc for node order (still safe to do globally)
    malloc1D(&node_order, NodeNum);

    eps0 = EpsT * 0.05; 
    eps1 = EpsT * 0.60; 
    eps2 = EpsT * 0.35; 

    p = 1.0 / (exp(eps1) + 1.0);
    V0 = p * (1.0 - p) / ((1.0 - 2.0 * p) * (1.0 - 2.0 * p));
    gamma = (1.0 - p) / (1.0 - 2.0 * p);
    lap_var_Y = (gamma * gamma) / (eps2 * eps2); 

    MakeRndPerm(node_order, NodeNum, NodeNum);

    int target_pairs = (PairNum != -1) ? std::min(PairNum, NodeNum / 2) : (NodeNum / 2);

    // Start OpenMP Parallel Region
    #pragma omp parallel reduction(+:cy4_num_ns)
    {
        // Initialize a thread-local Random Number Generator
        // We seed it with a base value + the thread ID so every thread has a unique sequence
        std::mt19937_64 local_rng(1776 + omp_get_thread_num());
        std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

        // Every thread gets its own private dense arrays
        double *local_a_mat_node1;
        double *local_a_mat_node2;
        malloc1D(&local_a_mat_node1, NodeNum);
        malloc1D(&local_a_mat_node2, NodeNum);

        for(int m = 0; m < NodeNum; m++) {
            local_a_mat_node1[m] = 0;
            local_a_mat_node2[m] = 0;
        }

        // Divide the loop iterations among the active threads
        #pragma omp for
        for(int pair_idx = 0; pair_idx < target_pairs; pair_idx++) {
            int k = pair_idx * 2;
            int node1 = node_order[k];
            int node2 = node_order[k+1];
            
            double cni2j, cnj2i, y_ij, x_ij, sigma2_ij, b_ij;
            double deg_i, deg_j;
            int true_edge, noisy_edge, v;
            double rnd;

            // Populate sparse map to thread-local dense array
            for (auto aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
                local_a_mat_node1[aitr->first] = 1;
            }
            for (auto aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
                local_a_mat_node2[aitr->first] = 1;
            }

            // Phase 2: Local One-sided CN from node1 to node2
            cni2j = 0.0;
            for (auto aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
                v = aitr->first;
                if (v == node2) continue; 
                
                true_edge = local_a_mat_node2[v]; 
                rnd = uniform_dist(local_rng); // Thread-safe random
                noisy_edge = 0;
                if (true_edge == 1 && rnd >= p) noisy_edge = 1;
                else if (true_edge == 0 && rnd < p) noisy_edge = 1;

                cni2j += ((double)noisy_edge - p) / (1.0 - 2.0 * p);
            }
            cni2j += LaplaceNoiseThreadSafe(gamma / eps2, local_rng); // Thread-safe noise

            // Phase 2: Local One-sided CN from node2 to node1
            cnj2i = 0.0;
            for (auto aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
                v = aitr->first;
                if (v == node1) continue; 
                
                true_edge = local_a_mat_node1[v]; 
                rnd = uniform_dist(local_rng); // Thread-safe random
                noisy_edge = 0;
                if (true_edge == 1 && rnd >= p) noisy_edge = 1;
                else if (true_edge == 0 && rnd < p) noisy_edge = 1;

                cnj2i += ((double)noisy_edge - p) / (1.0 - 2.0 * p);
            }
            cnj2i += LaplaceNoiseThreadSafe(gamma / eps2, local_rng); // Thread-safe noise

            // Phase 3: Dual-center average
            y_ij = (cni2j + cnj2i) / 2.0;

            true_edge = local_a_mat_node1[node2]; 
            rnd = uniform_dist(local_rng); // Thread-safe random
            noisy_edge = 0;
            if (true_edge == 1 && rnd >= p) noisy_edge = 1;
            else if (true_edge == 0 && rnd < p) noisy_edge = 1;
            x_ij = ((double)noisy_edge - p) / (1.0 - 2.0 * p);

            // Phase 0: Add Laplace noise to degrees
            deg_i = (double)a_mat[node1].size() + LaplaceNoiseThreadSafe(1.0 / eps0, local_rng);
            deg_j = (double)a_mat[node2].size() + LaplaceNoiseThreadSafe(1.0 / eps0, local_rng);
            
            // Variance estimation and 4-cycle debiasing
            sigma2_ij = ((deg_i + deg_j - 2.0 * x_ij) * V0 / 4.0) + lap_var_Y;
            b_ij = (y_ij * y_ij - y_ij - sigma2_ij) / 2.0;

            // Update local sum (OpenMP reduction handles thread summation)
            cy4_num_ns += b_ij;

            // FAST CLEANUP: Reset ONLY modified elements using thread-local arrays
            for (auto aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
                local_a_mat_node1[aitr->first] = 0;
            }
            for (auto aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
                local_a_mat_node2[aitr->first] = 0;
            }
        } // End of parallel for loop

        // Free thread-local memory before thread dies
        free1D(local_a_mat_node1);
        free1D(local_a_mat_node2);
    } // End of OpenMP Parallel Region

    free1D(node_order);
    
    // Scale up the sampled unbiased estimate to the whole graph
    double final_estimate = ((double)NodeNum * ((double)NodeNum - 1.0) / (4.0 * (double)target_pairs)) * cy4_num_ns;
    
    if (final_estimate < 0.0) {
        final_estimate = 0.0;
    }
    
    return final_estimate;
}

// ----------------------------------------------------------------------------
// SCN2P: Edge-LDP Sampled 2-Path (2-Star) Counting
// ----------------------------------------------------------------------------
double CalcSCN2P(map<int, int> *a_mat) {
    int *node_order;
    double eps1, eps2;
    double p, gamma;
    double st2_num_ns = 0.0; // Unbiased estimator sum for 2-stars

    // malloc for node order
    malloc1D(&node_order, NodeNum);

    // Empirical budget split for SCN2P (can be passed as parameters)
    eps1 = EpsT * 0.60;
    eps2 = EpsT * 0.40;

    p = 1.0 / (exp(eps1) + 1.0);
    gamma = (1.0 - p) / (1.0 - 2.0 * p);

    MakeRndPerm(node_order, NodeNum, NodeNum);

    int target_pairs = (PairNum != -1) ? std::min(PairNum, NodeNum / 2) : (NodeNum / 2);

    // Start OpenMP Parallel Region
    #pragma omp parallel reduction(+:st2_num_ns)
    {
        std::mt19937_64 local_rng(1776 + omp_get_thread_num());
        std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

        // Thread-private dense arrays for fast lookup
        double *local_a_mat_node1;
        double *local_a_mat_node2;
        malloc1D(&local_a_mat_node1, NodeNum);
        malloc1D(&local_a_mat_node2, NodeNum);

        for(int m = 0; m < NodeNum; m++) {
            local_a_mat_node1[m] = 0;
            local_a_mat_node2[m] = 0;
        }

        #pragma omp for
        for(int pair_idx = 0; pair_idx < target_pairs; pair_idx++) {
            int k = pair_idx * 2;
            int node1 = node_order[k];
            int node2 = node_order[k+1];
            
            double cni2j, cnj2i, y_ij;
            int true_edge, noisy_edge, v;
            double rnd;

            // Populate sparse map to thread-local dense array
            for (auto aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
                local_a_mat_node1[aitr->first] = 1;
            }
            for (auto aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
                local_a_mat_node2[aitr->first] = 1;
            }

            // Phase 2: Local One-sided CN from node1 to node2
            cni2j = 0.0;
            for (auto aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
                v = aitr->first;
                if (v == node2) continue; 
                
                true_edge = local_a_mat_node2[v]; 
                rnd = uniform_dist(local_rng);
                noisy_edge = 0;
                if (true_edge == 1 && rnd >= p) noisy_edge = 1;
                else if (true_edge == 0 && rnd < p) noisy_edge = 1;

                cni2j += ((double)noisy_edge - p) / (1.0 - 2.0 * p);
            }
            cni2j += LaplaceNoiseThreadSafe(gamma / eps2, local_rng);

            // Phase 2: Local One-sided CN from node2 to node1
            cnj2i = 0.0;
            for (auto aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
                v = aitr->first;
                if (v == node1) continue; 
                
                true_edge = local_a_mat_node1[v]; 
                rnd = uniform_dist(local_rng);
                noisy_edge = 0;
                if (true_edge == 1 && rnd >= p) noisy_edge = 1;
                else if (true_edge == 0 && rnd < p) noisy_edge = 1;

                cnj2i += ((double)noisy_edge - p) / (1.0 - 2.0 * p);
            }
            cnj2i += LaplaceNoiseThreadSafe(gamma / eps2, local_rng);

            // Phase 3: Dual-center average (First-order sum, NO quadratic correction)
            y_ij = (cni2j + cnj2i) / 2.0;
            st2_num_ns += y_ij;

            // Cleanup
            for (auto aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
                local_a_mat_node1[aitr->first] = 0;
            }
            for (auto aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
                local_a_mat_node2[aitr->first] = 0;
            }
        } // End of parallel for loop

        free1D(local_a_mat_node1);
        free1D(local_a_mat_node2);
    } // End of OpenMP Parallel Region

    free1D(node_order);
    
    // Scale up the sampled unbiased estimate to the whole graph
    // Scaling factor: n(n-1) / 2m
    double final_estimate = ((double)NodeNum * ((double)NodeNum - 1.0) / (2.0 * (double)target_pairs)) * st2_num_ns;
    
    if (final_estimate < 0.0) final_estimate = 0.0;
    
    return final_estimate;
}

// ----------------------------------------------------------------------------
// Imola SOTA Baseline: Laplace Degree Method for 2-Stars (USENIX Security 2021)
// ----------------------------------------------------------------------------
double CalcImolaSt(int* deg, double epsilon) {
    double est_total_2stars = 0.0;
    // Edge-LDP sensitivity for a node's degree is 1
    double b = 1.0 / epsilon;
    
    std::mt19937_64 local_rng(1776); // Single-threaded is fine since O(N) is extremely fast

    for (int i = 0; i < NodeNum; i++) {
        // 1. Local Perturbation: Add Laplace noise
        double noisy_deg = (double)deg[i] + LaplaceNoiseThreadSafe(b, local_rng);
        
        // 2. Server unbiased bias correction
        // E[d_i^2] = d_i^2 + 2b^2. We need to estimate d_i(d_i-1)/2 = (d_i^2 - d_i)/2
        double est_local_2stars = (noisy_deg * noisy_deg - 2.0 * b * b - noisy_deg) / 2.0;
        
        est_total_2stars += est_local_2stars;
    }

    if (est_total_2stars < 0.0) est_total_2stars = 0.0;
    return est_total_2stars;
}

// Calculate the clustering-coefficient
double CalcClstCoef(double tri_num_ns, double st2_num_ns){
	double clst_ns;

    if(tri_num_ns < 0) tri_num_ns = 0;
    if(st2_num_ns < 0) st2_num_ns = 0;
    if(st2_num_ns == 0) clst_ns = 1.0;
    else clst_ns = 3.0 * tri_num_ns / st2_num_ns;
    if(clst_ns > 1.0) clst_ns = 1.0;
    else if(clst_ns < 0.0) clst_ns = 0.0;

	return clst_ns;
}
int main(int argc, char *argv[])
{
    int all_node_num;
    int **node_order;
    map<int, int> *a_mat;           // adjacency matrix (for DP algorithms)
    vector<int> *adj_list;          // adjacency list (for fast true counts)
    map<int, int>::iterator aitr;
    map<int, int>::iterator aitr2;
    map<int, int>::iterator aitr3;
    int *deg;                                   // degree
    int *deg_lower;                             // degree in the lower-triangular part
    int max_deg;
    long long tri_num, st2_num;
    long long cy4_num;
    long long pa2_pow2;
    double clst;
    double tri_num_ns, sen_tri;
    double st2_num_ns, sen_st2;
    double clst_ns;
    double clst_re_ns, clst_l2_ns;
    double clst_re_ns_avg, clst_l2_ns_avg;
    double cy4_num_ns;
    double cy4_re_ns, cy4_l2_ns;
    double cy4_re_ns_avg, cy4_l2_ns_avg;
    int fix_perm;
    int itr;
    int i, j, k, l;
    string outdir;
    string outfile;
    char s[1025];
    char *tok;
    FILE *fp;
    clock_t start_clock, end_clock;
    double calc_time, calc_time_avg;
    double eclip_sum, tclip_sum;
    int eclip_num, tclip_num;

    // Initialization of Mersenne Twister
    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
    init_by_array(init, length);

    if (argc < 2) {
        printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon-delta (default: 1-8)] [#pairs (default:-1)] [#itr(-1) (default: 1)] [alg (default: 0)] ([bi (default: 0)]))]\n\n", argv[0]);
        printf("[EdgeFile]: Edge file\n");
        printf("[#nodes]: Number of nodes (-1: all)\n");
        printf("[epsilon-delta]: Privacy budgets (delta: -exponent)\n");
        printf("[#pairs]: Number of pairs (-1: all)\n");
        printf("[#itr(-1)]: Number of iterations (set #itr-1 to fix the permutation of nodes)\n");
        printf("[alg]: Algorithm (1: central triangle, 2(n): shuffle triangle (wedge) (n: numerical), 3(n)(-thr (default: 1)): shuffle triangle (wedge + ignore (thr * d_avg)) (n: numerical), 4(n)(-thr (default: 1)): shuffle triangle (wedge + ignore (thr * noisy d_avg) (n: numerical), 5: local triangle (wedge), 6(-smpl weight (default: 1)): local triangle (ARR), 7(-mu* (default: 1)): local triangle (2-rounds), 8(n): shuffle 4-cycle (n: numerical), 9: local 4-cycle)\n");
        printf("[bi]: Make a bipartite graph (1: yes, 0: no)\n");
        return -1;
    }

    EdgeFile = argv[1];

    NodeNum = -1;
    if (argc >= 3) NodeNum = atoi(argv[2]);

    EpsT = 1.0;
    Delta = 8.0;
    EpsT_s = "1";
    Delta_s = "8";
    if (argc >= 4){
        if((tok = strtok(argv[3], "-")) == NULL){
            printf("Error: incorrect [epsilon1-epsilon2]\n");
            exit(-1);
        }
        EpsT = atof(tok);
        EpsT_s = tok;

        if((tok = strtok(NULL, "-")) == NULL){
            printf("Error: incorrect [epsilon1-epsilon2]\n");
            exit(-1);
        }
        Delta = atof(tok);
        Delta_s = tok;
    }
    
    PairNum = -1;
    if (argc >= 5) PairNum = atoi(argv[4]);

    ItrNum = 1;
    fix_perm = 0;
    if (argc >= 6){
        tok  = strtok(argv[5], "-");
        ItrNum = atoi(tok);
        if((tok  = strtok(NULL, "-")) != NULL){
            if (strcmp(tok, "1") != 0){
                printf("Error: incorrect [#itr(-1)]\n");
                exit(-1);
            }
            else fix_perm = 1;
        }
    }

    Alg = 0;
    AlgPrm = 1;
    NumericalBound = 0;
    Alg_s = "0";
    if (argc >= 7){
        Alg_s = string(argv[6]);
        tok  = strtok(argv[6], "-");
        if(strlen(tok) == 1){
            Alg = atoi(tok);
        }
        else if(strlen(tok) == 2 && tok[1] == 'n'){
            Alg = tok[0] - '0';
            NumericalBound = 1;
        }
        else{
            printf("Error: incorrect [alg]\n");
            exit(-1);
        }

        if((tok  = strtok(NULL, "-")) != NULL){
            AlgPrm = atof(tok);
        }
    }

    Bip = 0;
    if (argc >= 8) Bip = atoi(argv[7]);

    // Total number of nodes --> all_node_num
    fp = FileOpen(EdgeFile, "r");
    for(i=0;i<2;i++) fgets(s, 1024, fp);
    all_node_num = atoi(s);
    fclose(fp);

    // malloc
    malloc2D(&node_order, ItrNum, all_node_num);

    // Use all nodes
    if (NodeNum == -1){
        NodeNum = all_node_num;
        for(j=0;j<NodeNum;j++) node_order[0][j] = j;
    }
    // Randomly generate the order of nodes --> node_order
    else{
        i = EdgeFile.find_last_of("/");
        outdir = EdgeFile.substr(0, i+1);
        outfile = outdir + "node-order_itr" + to_string(ItrNum) + ".csv";
        if(checkFileExistence(outfile)){
            fp = FileOpen(outfile, "r");
            for(j=0;j<all_node_num;j++){
                fgets(s, 1024, fp);
                strtok(s, ",");
                for(i=0;i<ItrNum;i++){
                    node_order[i][j] = atoi(strtok(NULL, ","));
                }
            }
            fclose(fp);
        }
        else{
            for(i=0;i<ItrNum;i++){
                MakeRndPerm(node_order[i], all_node_num, all_node_num);
            }
            fp = FileOpen(outfile, "w");
            for(j=0;j<all_node_num;j++){
                fprintf(fp, "%d,", j);
                for(i=0;i<ItrNum;i++) fprintf(fp, "%d,", node_order[i][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
        }

        // Use only the first permutation
        if (fix_perm){
            for(j=0;j<all_node_num;j++){
                for(i=1;i<ItrNum;i++) node_order[i][j] = node_order[0][j];
            }
        }
    }

    // Initialization
    malloc1D(&deg, NodeNum);
    malloc1D(&deg_lower, NodeNum);
    clst_re_ns_avg = clst_l2_ns_avg = 0.0;
    cy4_re_ns_avg = cy4_l2_ns_avg = 0.0;
    calc_time_avg = 0.0;

    // Output the header
    i = EdgeFile.find_last_of("/");
    outdir = EdgeFile.substr(0, i+1);
    if(Bip == 1){
        if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "b_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + "-1.csv";
        else outfile = outdir + "res_n" + to_string(NodeNum) + "b_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + ".csv";
    }
    else{
        if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + "-1.csv";
        else outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + ".csv";
    }
    fp = FileOpen(outfile, "w");
    if (abs(Alg - 5) == 4) {
        fprintf(fp, "#4cyc(true),#4cyc(est),#4cyc(rel-err),#4cyc(l2-loss),max_deg\n");
    } else {
        fprintf(fp, "#2hop(true),#2hop(est),#2hop(rel-err),#2hop(l2-loss),max_deg\n");
    }
    
    fclose(fp);

    // For each iteration
    for(itr=0;itr<ItrNum;itr++) {
        // Read edges for each iteration when NodeNum < all_node_num
        if(NodeNum < all_node_num || itr == 0){
            // Initialization
            a_mat = new map<int, int>[NodeNum];
            adj_list = new vector<int>[NodeNum]; // Initialize the fast adjacency list

            // Read edges from the edge file --> a_mat & adj_list
            ReadEdges(a_mat, adj_list, node_order[itr]);

            // Degree --> deg & max_deg
            max_deg = 0;
            st2_num = 0;
            for(i=0;i<NodeNum;i++){
                deg[i] = adj_list[i].size();
                if(max_deg < deg[i]) max_deg = deg[i];
                
                // #2-stars --> st2_num
                st2_num += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
            }

            // Degree --> deg_lower
            for(i=0;i<NodeNum;i++) deg_lower[i] = 0;
            for(i=0;i<NodeNum;i++){
                for (int neighbor : adj_list[i]){
                    if(neighbor < i) deg_lower[i] += 1;
                }
            }

            // Count 4-cycles in the original graph using the optimized vector structure
            if (fix_perm == 0 || itr == 0) {
                cy4_num = 0;
                pa2_pow2 = 0;
                
                // Fast O(1) counting using a global array and active tracking
                vector<int> pa2_cnt_fast(NodeNum, 0);
                vector<int> active_k; 
                active_k.reserve(1000); 
                
                for(i=0; i<NodeNum; i++){
                    for (int j : adj_list[i]) {
                        for (int k : adj_list[j]) {
                            // Exclude backward loops
                            if (k == i) continue; 
                            
                            // If this is the first time reaching k from i, mark it active
                            if (pa2_cnt_fast[k] == 0) {
                                active_k.push_back(k);
                            }
                            pa2_cnt_fast[k]++;
                        }
                    }
                    
                    // Iterate only over modified k's to calculate squared sum and quickly reset
                    for (int k : active_k) {
                        if (i < k) { // Prevent double counting (i < k)
                            long long count = pa2_cnt_fast[k];
                            pa2_pow2 += count * count;
                        }
                        pa2_cnt_fast[k] = 0; // O(1) fast reset
                    }
                    active_k.clear();
                }
                
                // We use the fact that (#2-paths)^2 = #2-stars + 4 * #4-cycles
                cy4_num = (pa2_pow2 - st2_num) / 4;
            }
        }

        /************************ Calculate sub-graph counts ************************/
        // Local 4-cycle
        // Calculate #4-cycles
        start_clock = clock();
        if (Alg == 9) {
            // Local 4-cycle (Imola RR)
            cy4_num_ns = CalImolaCy4(a_mat);
        } else if (Alg == 1) {
            // Proposed 4-cycle (CN4C)
            cy4_num_ns = CalCN4C(a_mat);
        } else if (Alg == 8) {
            // SOTA 2-Star (Imola Laplace)
            cy4_num_ns = CalcImolaSt(deg, EpsT); 
            cy4_num = st2_num; // Override the truth to 2-stars
        } else if (Alg == 2) {
            // Proposed 2-Star (SCN2P)
            cy4_num_ns = CalcSCN2P(a_mat);
            cy4_num = st2_num; // Override the truth to 2-stars
        } else if (Alg == 5) {
            // Proposed 4-cycle (OneRound)
            cy4_num_ns = CalOneR(a_mat);
            cy4_num = st2_num; // Override the truth to 2-stars
        }
        end_clock = clock();
        calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
        
        /**************************** Evaluate the loss *****************************/
        // relative error --> cy4_re_ns
        cy4_re_ns = fabs(cy4_num_ns - (double)cy4_num) / max((double)cy4_num, 0.001 * NodeNum);
        cy4_re_ns_avg += cy4_re_ns;
        // l2_loss --> cy4_l2_ns
        cy4_l2_ns = (cy4_num_ns - (double)cy4_num)*(cy4_num_ns - (double)cy4_num);
        cy4_l2_ns_avg += cy4_l2_ns;

        // Calculation time
        calc_time_avg += calc_time;

        /**************************** Output the results ****************************/
        fp = FileOpen(outfile, "a");
        if(MeasureTime){
            fprintf(fp, "%e,%e,%e,%e,%d,%e\n", 
            (double)cy4_num, cy4_num_ns, cy4_re_ns, cy4_l2_ns, max_deg, calc_time);
            printf("%e,%e,%e,%e,%d,%e\n", 
            (double)cy4_num, cy4_num_ns, cy4_re_ns, cy4_l2_ns, max_deg, calc_time);
        }
        else{
            fprintf(fp, "%e,%e,%e,%e,%d\n", 
            (double)cy4_num, cy4_num_ns, cy4_re_ns, cy4_l2_ns, max_deg);
        }
        fclose(fp);
        

        if(NodeNum < all_node_num || itr == ItrNum - 1){
            delete[] a_mat;
            delete[] adj_list; // Free the fast adjacency list to prevent memory leaks
        }
    }

    /************************* Output the results (AVG) *************************/
    
    cy4_re_ns_avg /= (double)ItrNum;
    cy4_l2_ns_avg /= (double)ItrNum;
    calc_time_avg /= (double)ItrNum;

    fp = FileOpen(outfile, "a");
    fprintf(fp, "function,AVG(rel-err),AVG(l2-loss)\n");
    fprintf(fp, "4-cycles,%e,%e\n", cy4_re_ns_avg, cy4_l2_ns_avg);
    if(MeasureTime) fprintf(fp, "CalcTime,%e\n", calc_time_avg);
    fclose(fp);
    

    // free all stuffs
    free2D(node_order, ItrNum);
    free1D(deg);
    free1D(deg_lower);

    return 0;
}