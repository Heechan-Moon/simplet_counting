#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<set>
#include<time.h>
#include<omp.h>

#include "file.h"
#include "CC.h"

using namespace std;

#include <sys/resource.h>

void printHeapMemoryUsage() {
    // Get heap memory usage using getrusage()
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    cout << usage.ru_maxrss << " KB" << endl;
}

int main(int argc, char **argv){
    //cout << "Default: ";printHeapMemoryUsage();
    int k=atoi(argv[1]); //k=4,5
    long long s=atoll(argv[2]); //s: sample size

    int tries=atoi(argv[3]); 
    string data=argv[4];
    int version=atoi(argv[5]);

    double duration;
    double start=clock();
    double start2=omp_get_wtime();

    
    ofstream myfile;
    myfile.open("../result/timestamp/CC_"+data+"_k"+to_string(k)+"_s"+to_string(s)+"_ver"+to_string(version)+".csv"); 

    File file(data);

    vector<vector<int>> Adj_list;
    
    file.simplicial_list();

    Adj_list=file.old_Adj; 
    //cout<<file.V<<endl;
    // file.print_graph(Adj_list, file.V); //
    //file.print_simplicies();
    
    myfile<<"reading,"<<(double)(clock()-start)/CLOCKS_PER_SEC<<"\n";
    //cout<<"  ===reading==="<<endl;
    start=clock();
    //cout << "Reading: ";printHeapMemoryUsage();
    
    for(int i=0; i<tries; i++){
        //myfile<<"try,"<<i<<"\n";
        //cout<<"try: "<<i<<endl;
        // make table
        start=clock(); start2 = omp_get_wtime(); 
        CC cc(k, file.V, file.simplices, file.vertex_to_simplices);

        cc.create_project_graph();
        if(version==1){
            cc.build(); 
        }else if(version==2){
            cc.build_without_memory();
        }

        
        //cout<<"building,"<<(double)(clock()-start)/CLOCKS_PER_SEC << " omptime " << (omp_get_wtime() - start2) << "\n";
        myfile<<"building,"<< (omp_get_wtime() - start2) << "\n";
        //cout<<"  ===building==="<<endl;

        //sampling
        //cout << "Building: ";printHeapMemoryUsage();
        start=clock(); start2=omp_get_wtime();
        if(version==1){
            cc.sampling(s);
        }else if(version==2){
            cc.sampling_without_memory(s);
        }
        //cout<<"sampling,"<<(double)(clock()-start)/CLOCKS_PER_SEC << " omptime " << (omp_get_wtime() - start2) << "\n";
        myfile<<"sampling,"<< (omp_get_wtime() - start2) << "\n";
        //cout<<"  ===sampling==="<<endl;
        //cout << "Sampling: ";printHeapMemoryUsage();
        //start=clock(); start2=omp_get_wtime();
        //cc.scan(file.simplices, file.vertex_to_simplices, s);
        
        //cout<<"scanning,"<<(double)(clock()-start)/CLOCKS_PER_SEC << " omptime " << (omp_get_wtime() - start2) << "\n";
        //myfile << "scaning,"<< (omp_get_wtime() - start2) << "\n";
        //cout<<"  ===scaning==="<<endl;

        //start=clock(); start2=omp_get_wtime();
        //cc.find_motifs(s, w);
        //cout<<"matching,"<<(double)(clock()-start)/CLOCKS_PER_SEC << " omptime " << (omp_get_wtime() - start2) << "\n";
        //myfile<<"matching,"<< (omp_get_wtime() - start2) << "\n";
        //cout<<"  ==matching=="<<endl;

        cc.normalize(s);

        //cc.print_final_result()
        cc.save_result(s, i, data, version);
        //cout << "Saving: ";printHeapMemoryUsage();
    } 
    myfile.close();
    //cout << "Final: ";printHeapMemoryUsage();
    return 0;
}