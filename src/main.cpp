#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<set>
#include<time.h>
#include<omp.h>
#include<math.h>
#include "file.h"
#include "RW.h"

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
    string data=argv[1];
    int k=atoi(argv[2]);
    long long s=atoll(argv[3]); //s: query size
    int tries=atoi(argv[4]); 
    
    int point = 20;

    double duration;
    double start=clock();
    double start2=omp_get_wtime();
    
    ofstream myfile;
    myfile.open("../result/timestamp/RW_"+data+"_k"+to_string(k)+"_s"+to_string(s)+"_.csv"); 

    File file(data);
    file.simplicial_list();
    
    //initialize
    myfile<<"reading,"<<(double)(clock()-start)/CLOCKS_PER_SEC<<"\n";
    start=clock();
    //cout << "Reading: ";printHeapMemoryUsage();
    for(int i=0; i<tries; i++){
        //myfile<<"try,"<<i<<"\n";
        start=clock(); start2 = omp_get_wtime(); 
        RW rw(k, file.V, file.simplices, file.vertex_to_simplices);
        
        rw.walk_uniform_online(s, point);
        //cout<<"sampling,"<<(double)(clock()-start)/CLOCKS_PER_SEC << " omptime " << (omp_get_wtime() - start2) << "\n";
        myfile<<"sampling,"<< (omp_get_wtime() - start2) << "\n";
        //cout<<"  ===sampling==="<<endl;
        //cout << "Sampling: ";printHeapMemoryUsage();
        rw.normalize();

        //cout<<"matching,"<<(double)(clock()-start)/CLOCKS_PER_SEC << " omptime " << (omp_get_wtime() - start2) << "\n";
        //myfile<<"matching,"<< (omp_get_wtime() - start2) << "\n";
        //cout<<"  ==matching=="<<endl;
        
        rw.save_result(s, i, data);
        //cout << "Saving: ";printHeapMemoryUsage();
    }
    myfile.close();
    //cout << "Final: ";printHeapMemoryUsage();
    return 0;
}