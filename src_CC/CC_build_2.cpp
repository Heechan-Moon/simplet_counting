#include "CC.h"
#include<vector>
#include<random>
#include<iostream>
#include<iomanip>
#include<map>
#include<bitset>
#include<math.h>
#include<typeinfo>
#include<omp.h>
#include<algorithm>

#define LIMIT 3372036854775807 // 9223372036854775807

using namespace std;

void CC::build_without_memory(){ //build color table

    //memoized_rng.resize(V);
    for(int i=0;i<V;i++){
        C[i].resize(cumul[k-1] + num[k-1]);
        //memoized_rng[i].resize(cumul[k-1] + num[k-1]);
        for(int sz=1;sz<=k;sz++){
            for(int j=cumul[sz-1];j<cumul[sz-1] + num[sz-1];j++){
                C[i][j].resize((int)bin_ones[sz-1].size());
                //memoized_rng[i][j].resize((int)bin_ones[sz-1].size());
            }
        }
        C[i][0][color[i]] = 1;
    }
    

    for(int curr_k=2;curr_k<=k;curr_k++){
        
        int num_subtasks = V * num[curr_k-1];
        //cout << "curr_k=" << curr_k << " start!" << endl;
        #pragma omp parallel for
        for(int task_i = 0; task_i < num_subtasks ;task_i++){
            
            int v = task_i / num[curr_k-1], dtype = cumul[curr_k-1] + (task_i % num[curr_k-1]); 
            if((curr_k < k) || (!color[v])){
                update_ver2(v, dtype);
            }
        }
        //cout << "curr_k=" << curr_k << " done." << endl;
    }
    rng_root.resize(V * num[k-1]);
    int decomp_from_2nd = cumul[k-1], decomp_to_2nd = cumul[k-1] + num[k-1];   
    count_type cumul_counts = 0;
    double _start2 = omp_get_wtime();
    for(int v=0;v<V;v++){
        for(int i=0;i<num[k-1];i++){
            cumul_counts += C[v][cumul[k-1] + i][0];
            rng_root[v * num[k-1] + i] = cumul_counts;
        }
    }
    //cout << (omp_get_wtime() - _start2) << " New:" << cumul_counts << endl;
}


void CC::update_ver2(int v, int treelet_id){ //update color table
    for(auto &u: adj_list[v]){
        for(auto &decomp_case: D[treelet_id]){ //candidate decomposition of T as T1 and T2
            int main_num_nodes = D[decomp_case.main_id][0].num_nodes;
            for(auto &cs: color_pairs[main_num_nodes][decomp_case.num_nodes - main_num_nodes]){
                C[v][treelet_id][color_idx[cs.first | cs.second]] += C[v][decomp_case.main_id][color_idx[cs.first]] * C[u][decomp_case.sub_id][color_idx[cs.second]];
                //memoized_rng[v][treelet_id][color_idx[cs.first | cs.second]].capacity += 1;
            }
        }
    }
    for(auto &target_color: bin_ones[D[treelet_id][0].num_nodes-1]){
        C[v][treelet_id][color_idx[target_color]] /= D[treelet_id][0].deg_root;
    }
}

