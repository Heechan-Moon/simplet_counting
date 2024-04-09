#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<fstream>
#include<iostream>
#include"file.h"
#include<algorithm>

using namespace std;

File::File(string d){

    filepath="../dataset/";
    //proj="-proj-graph";

    filename_nv=filepath+d+"/"+d+"-nverts.txt";
    filename_sim=filepath+d+"/"+d+"-simplices.txt";
    //filename_proj=filepath+d+proj+"/"+d+proj+".txt";
   
}

void File::simplicial_list(){ 
    
    get_simplices_list();
    get_adjacency_list();

    return;
}

void File::get_simplices_list(){
    ifs.open(filename_nv); ifs2.open(filename_sim);
    
    if(ifs.is_open() && ifs2.is_open()){
        int cnt = 0;
        while(getline(ifs, s)){
            num=stoi(s);
            vector<int> temp;
            for(int i=0; i<num; i++){
                getline(ifs2, s2);
                int num2=stoi(s2) - 1;
                if(V <= num2) V = num2 + 1;
                temp.push_back(num2);
            }
            if(temp.size()>1){
                sort(temp.begin(), temp.end());
                temp.erase(unique(temp.begin(),temp.end()),temp.end());
                simplices.push_back(temp);
            }
        }
    }
    ifs.close();  ifs2.close();
    return;
}

void File::get_adjacency_list(){
    vertex_to_simplices.resize(V);
 
    int size=simplices.size();
    for(int i=0; i<size; i++){
        int l=simplices[i].size();
        for(int j=0; j<l; j++){
            vertex_to_simplices[simplices[i][j]].push_back(i);
        }
        //projected_graph(i);
    }

    return;
}


void File::print_simplices(){
    for(int i=0; i<simplices.size(); i++){
        for(int j=0; j<simplices[i].size(); j++){
            cout<<simplices[i][j]<<" ";
        }
        cout<<endl;
    }
}


void File::addEdge(vector<vector<int>> *adj, int u, int v){
    (*adj)[u].push_back(v);
    (*adj)[v].push_back(u);
}


void File::print_graph(vector<vector<int>> adj, int V){
    for (int v =0; v < V; ++v) {
        cout << "vertex " << v+1 << " ";
        for (auto x : adj[v])
            cout << "-> " << x;
        cout<<endl;
    }
}

/*
void File::weighted_edge_list(){
    tot_Adj.resize(V);
    ifs3.open(filename_proj);
    if(ifs3.is_open()){
        while(getline(ifs3, s)){
            istringstream ss(s);
            vector<int> temp;
            while(getline(ss, s, ' ')){                
                temp.push_back(stoi(s));        
            } 
            edges.push_back(temp);
            addEdge(&tot_Adj, temp[0], temp[1]);
            if(V<temp[1]){V=temp[1];}
        }
    }
    ifs3.close();  
    return;
}*/