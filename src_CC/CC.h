#ifndef CC_H
#define CC_H
#include<string>
#include<vector>
#include<map>
#include<tuple>
#include<set>
#include<math.h>
#include<unordered_map>
#include <atomic>
#include <queue>

using namespace std;
using count_type = __int128_t;
//using count_type = long long;

class CC {
    struct decomp_type{
        int num_nodes, main_id, sub_id, deg_root;
    };
    struct queue_type{
        int v, treelet_id, colors;
    };
    struct rng_type{
        struct info_type{
            int main_tid, main_colors, sub_v, sub_tid, sub_colors;
        };
        int capacity;
        vector<info_type> infos;
        vector<count_type> cumcnts;
    };
    int n_threads;
    int k; 
    int V; 
    vector<vector<int>> simplices;
    vector<vector<int>> vertex_to_simplices;
    vector<vector<int>> adj_list;
    
    void initialize_motifs();

    count_type total_treelet=0;
    
    vector< vector<decomp_type> > D;
    
    int sum=0;
    int colors;
    int weight;    

    int num[6]={1, 1, 2, 4, 9, 20}; //number of rooted trees of n nodes
    int diff[6]={0, 0, 0, 1, 4, 12}; // cumul-i
    int cumul[6]={0,1,2,4,8,17}; //for tree 17-> (4, 9)
    
    vector<int> num_ones;
    vector< vector<int> > bin_ones;
    vector<int> color_idx;

    vector< vector< vector<count_type> > > C;
    vector< vector< vector<rng_type> > > memoized_rng;

    
    //vector< vector<queue_type> > ques;

    void create_tree_table();
    void update(int v, int treelet_id); 
    void update_ver2(int v, int treelet_id); 
    vector< vector< vector< pair<int, int> > > > color_pairs;
    vector<int> color;

    //sample
    vector<count_type> rng_root;

    //void iterative_sample(int t);    
    void make_sample_and_update(int tid, long long idx);
    void make_sample_and_update_ver2(int tid, long long idx);
    vector<int> scan(vector<int> sample, int tid, long long sample_idx);
    //vector< vector<int> > sample;
    vector< vector<int> > seen_vertices, intersect_counts;
    vector< vector<long long> > seen_simplices;
    //scan
    //vector< vector< pair<int, int> > > maximal_simplices;
    
    //motif
    int factorial(int n);
    vector<int> num_spanning; //vector<double> density;

    vector<vector<int>> permuted_code;
    int tot_motif=0;
    int n_simplets=0;
    
    unordered_map<long long, int> hash;
    vector<int> spanning_cnts;

    count_type normalize_constant=0;

    
    void load_spannings();
    void load_motifs();
    //void final_motifs(int s);
    
    public:
        CC(int K, int V, vector<vector<int>> simplices, vector<vector<int>> v_to_sim);
        void build();
        void build_without_memory();
        void sampling(long long s);
        void sampling_without_memory(long long s);
        void create_project_graph();
        //void find_motifs(int s, int w);
        //void print_gmotifs();
        //void print_hashkey();
        //void print_final_result();
        void save_result(long long s, int trial, string task, int version);
        void normalize(long long s);
        vector<double> simplets;

};

#endif    