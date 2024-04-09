#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<fstream>
#include<iostream>
#include"RW.h"
#include<algorithm>
#include<random>
#include<utility>
#include<queue>
#include<deque>
#include<set>
#include<unordered_set>
#include<omp.h>
#include<iterator>


RW::RW(int k, int V, vector<vector<int>> simplices, vector<vector<int>> v_to_sim){
    
    this -> V = V;
    this -> k = k;
    this -> simplices = simplices;
    this -> vertex_to_simplices = v_to_sim;
    
    this -> n_threads = omp_get_max_threads();

    initialize_motifs();

    states.resize(n_threads);
    trans_prob.resize(n_threads);
    prior.resize(n_threads);
    X.resize(n_threads);
    degree_temp.resize(n_threads);
    neighbor_temp.resize(n_threads);

    int _pow = (1<<k);
    num_ones.resize(_pow);
    bin_ones.resize(k);
    for(int i=1; i<_pow; i++){
        num_ones[i] = num_ones[i & (i-1)] + 1;
        //printf(" %d", num_ones[i]);
        bin_ones[num_ones[i]-1].push_back(i);
    }

    
    initialize_states();

    return;
}


void RW::initialize_states(){

    Sigma = simplices.size();

    random_device rd;
    mt19937 gen(rd());
    int num_sim = simplices.size();
    int e;
    int size;

    //#pragma omp parallel for
    for(int tid=0; tid<n_threads; tid++){

        uniform_int_distribution<int> dis(0, num_sim-1);   

        while(true){
            e = dis(gen);
            size = simplices[e].size(); 
            if(size > 1) break;
        }

        uniform_int_distribution<int> dis2(0, size-1);
        uniform_int_distribution<int> dis3(0, size-2);

        int a = dis2(gen);
        int b = dis3(gen);

        if(b >= a) b++;
        
        int st = simplices[e][a];
        int end = simplices[e][b];
        
        states[tid].push(make_pair(st, end));

        double C = 0;
        vector<int> cand;

        unordered_set<int> nei;
        for(auto &sim: vertex_to_simplices[st]){
            for(auto &u: simplices[sim]){if(u!=st) nei.insert(u);}
        }
        vector<int> nei_a(nei.begin(), nei.end());
        int temp_size = nei_a.size();
        nei.clear();
        for(auto &sim: vertex_to_simplices[end]){
            for(auto &u: simplices[sim]){if(u!=end) nei.insert(u);}
        }
        temp_size += nei.size();
        vector<int> nei_b(nei.begin(), nei.end());
        degree_temp[tid].push(temp_size-2);
        neighbor_temp[tid] = make_pair(nei_a, nei_b);

    }
   
    return;
}
