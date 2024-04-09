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

void CC::sampling_without_memory(long long s){
    long long tid_s = 0;
    int remains = s % n_threads;
   
    #pragma omp parallel for
    for(int tid=0; tid<n_threads; tid++){
        
        if(tid == n_threads-1) tid_s = s / n_threads + s % n_threads;
        else tid_s = s / n_threads;
        
        for(long long i=0; i<tid_s; i++){
            make_sample_and_update_ver2(tid, i);
        }
    }

    return;
}



void CC::make_sample_and_update_ver2(int tid, long long idx){

    thread_local random_device rd;   
    thread_local mt19937 gen(rd());
    thread_local uniform_real_distribution<double> dis(0, 1);
    rng_type non_memoized_rng;
    vector<int> sample;
    vector<queue_type> ques;
    sample.resize(k);
    ques.resize(k * 2 + 1);
    
    count_type where = (count_type)(dis(gen) * rng_root[(int)rng_root.size() - 1]) + 1;
    int lo = 0, hi = (int)rng_root.size() - 1, mid, sampled_idx = -1;

    while(lo <= hi){
        mid = (lo + hi) >> 1;
        if(where <= rng_root[mid]){
            sampled_idx = mid; hi = mid  - 1;
        }else{
            lo = mid + 1;
        }
    }
    int qhead = 0, qtail = 0;
    ques[qtail] = {sampled_idx / num[k-1], cumul[k-1] + (sampled_idx % num[k-1]), (1 << k) - 1}; qtail++;
    int pos = 0;

    while(qhead < qtail){
        queue_type root = ques[qhead]; qhead++;
        // cout << "curr top: " << root.v << " " << root.treelet_id << " " << C_new[root.v][root.treelet_id][root.colors] << endl;
        if(C[root.v][root.treelet_id][color_idx[root.colors]] < 1){
            break;
        }
        //cout << D[root.treelet_id][0].num_nodes << endl;
        if(D[root.treelet_id][0].num_nodes == 1){
            //cout << "curr top: " << root.v << " " << root.treelet_id << " " << root.colors << " " << C_new[root.v][root.treelet_id][root.colors] << endl;
            sample[pos] = root.v;
            pos++;
            continue;
        }


        //while (lock.test_and_set(memory_order_acquire));
        // allocate memory
    
        //lock.clear(memory_order_release);

        count_type curr_cnts = 0;
        //int infos_cnt = 0;
        for(auto &u: adj_list[root.v]){
            for(auto &decomp_case: D[root.treelet_id]){ //candidate decomposition of T as T1 and T2
                int main_num_nodes = D[decomp_case.main_id][0].num_nodes;
                for(auto &cs: color_pairs[main_num_nodes][decomp_case.num_nodes - main_num_nodes]){
                    if(root.colors ^ (cs.first | cs.second)) continue;
                    curr_cnts += C[root.v][decomp_case.main_id][color_idx[cs.first]] * C[u][decomp_case.sub_id][color_idx[cs.second]];

                    non_memoized_rng.infos.push_back({decomp_case.main_id, cs.first, u, decomp_case.sub_id, cs.second});
                    non_memoized_rng.cumcnts.push_back(curr_cnts);
                    //infos_cnt += 1;
                
                    
                }
            }
        }
        
        int choice_capacity = non_memoized_rng.infos.size();

        count_type where = (count_type)(dis(gen) * non_memoized_rng.cumcnts.back()) + 1;
        
        int lo = 0, hi = choice_capacity - 1, mid, sampled_idx = -1;
        while(lo <= hi){
            mid = (lo + hi) >> 1;
            if(where <= non_memoized_rng.cumcnts[mid]){
                sampled_idx = mid; hi = mid  - 1;
            }else{
                lo = mid + 1;
            }
        }

        //int sampled_idx = bsearch(memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].cumcnts, where);
        rng_type::info_type chosen_info = non_memoized_rng.infos[sampled_idx];
        ques[qtail] = {root.v, chosen_info.main_tid, chosen_info.main_colors}; qtail++;
        ques[qtail] = {chosen_info.sub_v, chosen_info.sub_tid, chosen_info.sub_colors}; qtail++;

        non_memoized_rng.infos.clear();
        non_memoized_rng.cumcnts.clear();
    }
    
    vector<int> maximal_simplices = scan(sample, tid, idx);

    long long hashed_code = 0;
    for(auto &key: maximal_simplices){
        hashed_code += (1LL << (key - 1));
    }

    int answer = hash[hashed_code];
  
    simplets[answer]+= 1;  
    return;

}
