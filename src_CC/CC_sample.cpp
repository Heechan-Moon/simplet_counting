#include "CC.h"
#include<vector>
#include<random>
#include<iostream>
#include<map>
#include<math.h>
#include <sstream>
#include<algorithm>
#include <queue>
#include <omp.h>
#include <mutex>

using namespace std;
int bsearch(const vector<count_type> &target_vec, const count_type &where){
    int lo = 0, hi = (int)target_vec.size() - 1, mid, target_idx = -1;
    while(lo <= hi){
        mid = (lo + hi) >> 1;
        //cout << lo << " " << hi << " " << target_vec[mid] << " " << where << endl;
        if(where <= target_vec[mid]){
            target_idx = mid; hi = mid  - 1;
        }else{
            lo = mid + 1;
        }
    }
    return target_idx;
}

void CC::make_sample_and_update(int tid, long long idx){

    thread_local random_device rd;   
    thread_local mt19937 gen(rd());
    thread_local uniform_real_distribution<double> dis(0, 1);
    
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
        int choice_capacity = memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].capacity;
        
        #pragma omp critical (update_rng)
        {  
            if(memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].infos.size() == 0){
                
                //while (lock.test_and_set(memory_order_acquire));
                // allocate memory
                
                memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].infos.resize(choice_capacity);
                memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].cumcnts.resize(choice_capacity);
                //lock.clear(memory_order_release);

                count_type curr_cnts = 0;
                int infos_cnt = 0;
                for(auto &u: adj_list[root.v]){
                    for(auto &decomp_case: D[root.treelet_id]){ //candidate decomposition of T as T1 and T2
                        int main_num_nodes = D[decomp_case.main_id][0].num_nodes;
                        for(auto &cs: color_pairs[main_num_nodes][decomp_case.num_nodes - main_num_nodes]){
                            if(root.colors ^ (cs.first | cs.second)) continue;
                            curr_cnts += C[root.v][decomp_case.main_id][color_idx[cs.first]] * C[u][decomp_case.sub_id][color_idx[cs.second]];

                            memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].infos[infos_cnt] = {decomp_case.main_id, cs.first, u, decomp_case.sub_id, cs.second};
                            memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].cumcnts[infos_cnt] = curr_cnts;
                            infos_cnt += 1;
                            if(infos_cnt > choice_capacity){
                                cout << infos_cnt << " " << choice_capacity << endl;
                            }
                            
                        }
                    }
                }
                if(choice_capacity != infos_cnt){
                    cout << "CAPACITY??" << endl;
                }
            }
        }

        count_type where = (count_type)(dis(gen) * memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].cumcnts.back()) + 1;
        
        int lo = 0, hi = choice_capacity - 1, mid, sampled_idx = -1;
        while(lo <= hi){
            mid = (lo + hi) >> 1;
            if(where <= memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].cumcnts[mid]){
                sampled_idx = mid; hi = mid  - 1;
            }else{
                lo = mid + 1;
            }
        }

        //int sampled_idx = bsearch(memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].cumcnts, where);
        rng_type::info_type chosen_info = memoized_rng[root.v][root.treelet_id][color_idx[root.colors]].infos[sampled_idx];
        ques[qtail] = {root.v, chosen_info.main_tid, chosen_info.main_colors}; qtail++;
        ques[qtail] = {chosen_info.sub_v, chosen_info.sub_tid, chosen_info.sub_colors}; qtail++;
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

void CC::sampling(long long s){
    long long tid_s = 0;
    int remains = s % n_threads;
   
    #pragma omp parallel for
    for(int tid=0; tid<n_threads; tid++){
        
        if(tid == n_threads-1) tid_s = s / n_threads + s % n_threads;
        else tid_s = s / n_threads;
        
        for(long long i=0; i<tid_s; i++){
            make_sample_and_update(tid, i);
        }
    }

    return;
}
