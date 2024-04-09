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
#include<map>

void RW::walk_uniform_online(long long s, int point){
    #pragma omp parallel for
    for(int tid=0; tid<n_threads; tid++){
        for(int i=0; i<k-2; i++){ 
            next_uniform_online(tid);
        }
    }
    long long tid_s = s / n_threads;
    #pragma omp parallel for
    for(int tid=0; tid<n_threads; tid++){
        for(long long i=0; i<tid_s; i++){
            //printf("%d-%d\n", i, tid);
            //print_states(states[tid]);
            queue<pair<int, int>> states_copy = states[tid];
            unordered_set<int> valid_sample;
            while (!states_copy.empty()) {
                pair<int, int> element = states_copy.front();
                states_copy.pop();
                valid_sample.insert(element.first);
                valid_sample.insert(element.second);
            }
            if(valid_sample.size()==k){update_simplet_counts(tid, i);}
            
            next_uniform_online(tid);   
            states[tid].pop();
            
            degree_temp[tid].pop();//add
        }
    }

    return;

}


void RW::next_uniform_online(int tid){ // next state (edge)

    random_device rd;
    mt19937 gen(rd());
    pair<int, int> curr = states[tid].back();
    int a = curr.first; // vertex
    int b = curr.second; // vertex
    int st; int old; int end;

    //to here
    uniform_int_distribution<int> dis(0, degree_temp[tid].back()-1);
    int temp = dis(gen);
    if(temp < neighbor_temp[tid].first.size()-1){
        st = a; old = b;
    }
    else{ 
        st = b; old = a;
        swap(neighbor_temp[tid].first, neighbor_temp[tid].second);
    }
    //printf("\n===========\n");
    //printf("start: %d\n", st);
    //3. convert unordered_set into vector
    uniform_int_distribution<int> dis2(0, neighbor_temp[tid].first.size() - 1);
    
    // choose a random node
    //for(auto &v: proj_graph[st]){
     //   printf("%d-", v);
   // }printf("\n");
    while(true){
        temp = dis2(gen);
        end = neighbor_temp[tid].first[temp]; 
        if(old!=end) break;
    }

    //change from here
    unordered_set<int> nei;
    for(auto &sim: vertex_to_simplices[end]){
        for(auto &u: simplices[sim]){if(u!=end) nei.insert(u);}
    }
    vector<int> nei_new(nei.begin(), nei.end());

    states[tid].push(make_pair(st, end));
    degree_temp[tid].push(neighbor_temp[tid].first.size()+nei_new.size()-2);
    neighbor_temp[tid] = make_pair(neighbor_temp[tid].first, nei_new);

    return;
}


void RW::update_simplet_counts(int tid, long long sample_idx){
    vector<int> ss;
    double amount = 1;
    queue<int> degree_copy = degree_temp[tid]; 
    queue<pair<int, int>> q_copy = states[tid];
    vector<pair<int,int>> q_copy_vector;
    
    queue<pair<int, int>> q_copy_temp = q_copy;
    

    // Printing all elements in the copy of the queue
    while (!q_copy_temp.empty() ) {
        pair<int, int> element = q_copy_temp.front();
        q_copy_vector.push_back(element);
        q_copy_temp.pop();
    }
    int idx = 0;
    while (!q_copy.empty())
    {
        pair<int, int> edge = q_copy.front();
        if(idx){
            ss.push_back(edge.second);
        }else{
            ss.push_back(edge.first);
            ss.push_back(edge.second);
        }
        q_copy.pop();
        idx ++;
    }    

    sort(ss.begin(), ss.end());
    ss.erase(unique(ss.begin(),ss.end()),ss.end());

    degree_copy.pop();
    while(degree_copy.size()!=1){
        amount *= (double)(1.0/degree_copy.front());
        degree_copy.pop();
    }

    vector<int> maximal_simplices = scan(ss, tid, sample_idx);

    long long hashed_code = 0;
    for(auto &key: maximal_simplices){
        hashed_code += (1LL << (key - 1));
    }

    int answer = hash[hashed_code];

    simplets[answer]+= 1/(amount);  
    s++;
    return;
}

void RW::print_states(queue<pair<int, int>> copy){

    cout<<"  states: ";
    
    while (!copy.empty())
    {
        pair<int, int> temp = copy.front();
        printf("(%d-%d) ", temp.first, temp.second);
        copy.pop();
    }
    printf("\n");
}

vector<int> RW::scan(vector<int> sample, int tid, long long i){

    vector<int> maximal_simplices;
    vector<int> intersect_counts;
    intersect_counts.resize((1 << k));
    
    //cout << "n_nodes: " << ((int)(vertex_to_simplices.size())) << ", n_simplices: " << (int)simplices.size() << endl;

    for(int j=0;j<k;j++){
        seen_vertices[tid][sample[j]] = (1 << j);
    }
    for(auto &v: sample){
        for (auto &sid: vertex_to_simplices[v]) {
            if(seen_simplices[tid][sid] == i+1) continue;
            seen_simplices[tid][sid] = i+1;
            int intersect_state = 0;
            for(auto &u: simplices[sid]){
                intersect_state |= seen_vertices[tid][u];
            }
            if(num_ones[intersect_state] >= 2){
                intersect_counts[intersect_state] += 1;
            }
        }
    }
    for(int j=0;j<k;j++){
        seen_vertices[tid][sample[j]] = 0;
    }

    for(int j=(1 << k)-1;j>=0;j--){
        if(intersect_counts[j] > 0){
            maximal_simplices.push_back(j);
            //cout<<j<<" ";
            for(int jj = 0; jj < j; jj++){
                if((j & jj) == jj) intersect_counts[jj] = 0;
            }
        }
    }


    //printf("==========\n"); for(auto &v: sample){ printf("%d ", v); } printf("\n"); for(auto &v: maximal_simplices){ printf("%d ", v); } printf("\n");
    return maximal_simplices;
}

