#ifndef RW_H
#define RW_H

#include<utility>
#include<queue>
#include<deque>
#include<unordered_map>
#include <string>

using namespace std;

class RW{    
    
    int k;
    int V;
    int Sigma;

    vector<vector<int>> simplices;
    vector<vector<int>> vertex_to_simplices;
    vector<int> induced_subgraph;
    vector<queue<pair<int, int>>> states;
    vector<queue<double>> n2v_weights;
    vector<int> weight_table;
    
    vector<queue<double>> trans_prob;
    
    //for uniform sampling
    vector<vector<int>> proj_graph;
    vector<int> degree;
    vector<int> cumul_degree;

    int n_threads;
    long long s;
    void initialize_states();    
    void initialize_motifs();
    void load_coeffs();
    void load_motifs();

    void next_uniform_online(int tid);
    void update_simplet_counts(int tid, long long sample_idx);
    
    vector<int> scan(vector<int> sample, int tid, long long sample_idx);

    //uniform_online
    vector<queue<int>> degree_temp;
    vector<pair<vector<int>, vector<int>>> neighbor_temp;

    //add below
    vector<int> num_ones;
    vector< vector<int> > bin_ones;

    vector<double> stationary_dist;
    
    vector<pair<vector<int>, int>> prior;
    vector<queue<double>> X;
    
    //motif
    vector< vector<int> > seen_vertices, seen_simplices, intersect_counts;
    
    int factorial(int n);
    vector<int> num_coeff; 

    vector<vector<int>> permuted_code;
    int tot_motif=0;
    int n_simplets=0;
    
    unordered_map<long long, int> hash;
    vector<int> coeff_cnts;

    public:
        void print_states(queue<pair<int, int>> copy);

        RW(int k, int V,  vector<vector<int>> simplices, vector<vector<int>> v_to_sim);
        

        void normalize();
        void walk_uniform_online(long long s, int point);
        
        vector<double> simplets;
        void save_result(long long s, int trial, string task);

};


#endif

