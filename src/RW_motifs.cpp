#include "RW.h"
#include<vector>
#include<iostream>
#include<map>
#include<algorithm>
#include<math.h>
#include<fstream>
#include<string>
#include<iterator>
#include<sstream>
#include<set>
#include<cmath>
#include <unordered_set>

void RW::initialize_motifs(){
    vector<int> base;
    for(int i=0; i<k; i++) base.push_back(i);
    int perm_cnt = 0; 
    do{
        permuted_code.push_back(vector<int>(1 << k, 0));
        for(int j=0;j<(1 << k);j++){
            for(int digit=0; digit<k; digit++){
                if(j & (1 << digit)) permuted_code[perm_cnt][j] += (1 << base[digit]);
            }
        }
        perm_cnt += 1;
    }while (next_permutation(base.begin(), base.end()));
    
    load_coeffs();
    load_motifs();
    simplets.resize(n_simplets, 0);
    
    seen_vertices.resize(n_threads); seen_simplices.resize(n_threads);// intersect_counts.resize(s);
    for(int i=0;i<n_threads;i++){
        seen_vertices[i].resize((int)vertex_to_simplices.size());
        seen_simplices[i].resize((int)simplices.size());
    }

    return;
}

void RW::normalize(){

    for(int i=0; i<n_simplets; i++){
        simplets[i] = (simplets[i] / (double)coeff_cnts[i]) / (double)s;
        
    }

    return;
}


void RW::load_coeffs(){
    num_coeff.clear();
    
    ifstream ifs; string s;
    ifs.open("../motif/rw_coeff_"+to_string(k)+".txt");
    if(ifs.is_open()){ 
        while(getline(ifs, s)){
            num_coeff.push_back(stoi(s));
        }
    }
    ifs.close(); 
    
    return;
}



void RW::load_motifs(){

    ifstream ifsm; string sm;
    ifstream ifse; string se;
    ifstream ifsk; string sk;
    ifsm.open("../motif/nummotifs_"+to_string(k)+".txt");
    ifse.open("../motif/numedges_"+to_string(k)+".txt");
    ifsk.open("../motif/key_"+to_string(k)+".txt");
    string s; 
    int zero=int('0');
    int perm_cnt = (int)permuted_code.size();
    n_simplets = 0;
    int gmotif_idx = 0;

    coeff_cnts.clear();
    vector<int> closed;
    //cout << "perm_cnt: " << perm_cnt << endl;
    if(ifsm.is_open() && ifse.is_open() && ifsk.is_open()){ // && ifsv.is_open()
        while(getline(ifsm, sm)){ // number of motifs
            int sm_size=stoi(sm); 
            for(int i=0; i<sm_size; i++){ //for each graph-motif 
                if(getline(ifse, se)){ //number of edges in each motifs
                    int se_size=stoi(se);
                    closed.clear();
                    for(int j=0; j<se_size; j++){
                        if( getline(ifsk, sk)){ //getline(ifsv, sv) &&
                            int raw_key = stoi(sk); //, raw_value = stoi(sv);
                            //if(raw_value > 0){
                                int simplex_code = 0;
                                while(raw_key){
                                    simplex_code += (1 << ((raw_key % 10) - 1));
                                    raw_key /= 10;
                                }
                                closed.push_back(simplex_code);
                            //}
                        }
                    }
                    for(int i=0;i<perm_cnt;i++){
                        long long hashed = 0;
                        for(auto &c_simplex: closed){
                            hashed |= (1LL << (permuted_code[i][c_simplex] - 1));
                        }
                        hash[hashed] = n_simplets;
                    }
                    coeff_cnts.push_back(num_coeff[gmotif_idx]);
                }
                n_simplets += 1;
            }
            gmotif_idx += 1;
        }
    }

    ifsm.close();
    ifse.close();
    ifsk.close();
    return;
}


int RW::factorial(int n){
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;

}

