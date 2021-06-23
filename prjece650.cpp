#include <minisat/core/Solver.h>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <list>
#include <pthread.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <map>
#include <set>
#include <numeric>
#include <cmath>
#include <chrono>

using namespace std;
//vector<int> adj1[100];
std::vector<int>temp1dedges;
std::vector<int>vertexcoverfinalresult;
//int vcresult[25];
std::vector<int> approxvc2_result;
std::vector <int> approxvc1_result;
int edgepair=0;
int Vertices=0;
string TempEdges;
int statusflag=0;
string inputval;
std::vector<int>vertexcoverresult;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

int eof_all=0;
int timeout=0;
int result_flag=0;

struct timespec time_s1, time_e1,time_s2, time_e2,time_s3, time_e3;
clockid_t threadClockId1, threadClockId2, threadClockId3;
double cnf_time, vc1_time, vc2_time;

class Param_Class_{
	public:
		std::vector<int> get_vertexCover(vector<vector<int>> adj_list, int V) ;
		void get_input(string line, int &V);
		std::vector<int> vcres;
		void print();

		void APPROX_VC_1(int V, unsigned int no_edges);
		std::vector<int> apvc1res;

		int approx_vc_2_recurrsivefunc(std::vector<int>approxvc2_list, int edgepair);
		int APPROX_VC_2(int edgepair);
		std::vector<int> apvc2res;

};

struct _threadData_{
    Param_Class_ inputgraph;
};

std::vector<std::string> splitfunction(std::string str, char splitby)
{
	std::stringstream ss(str);
	std::string item;
	std::vector<std::string> splittedStrings;
	while (std::getline(ss, item, splitby))
		{
			splittedStrings.push_back(item);
		}
	return splittedStrings;

}

void extractIntegersFromWords(string str)
{
	stringstream functionstr;
	functionstr << str;
	string tempstr;
	int found;

	while (!functionstr.eof()) {
		functionstr >> tempstr;

		if (stringstream(tempstr) >> found)
		{
			temp1dedges.push_back(found);
		}
		// tempstr = "";
	}


}

// --------------Logic for approx-vc-1-----------------------//
void Param_Class_::APPROX_VC_1(int V,unsigned int no_edges){
	std::vector<vector<int>> edges;
	for (unsigned int i=0; i<temp1dedges.size()-1; i++){
            if (i%2==0){
                //cout<<"i="<<i<<" "<<temp1dedges[i]<<" "<<temp1dedges[i+1]<<endl;
                edges.push_back({temp1dedges[i],temp1dedges[i+1]});
            }
	}

    int v1,v2;
    int vertex=-1;
    bool all;

    do{
        vector <int> adj_list[V];
        vector <int> nodes;
        int neighbours[V][2];

        int maximum=0;

        for (unsigned int i=0;i<no_edges;i++){
            if(edges[i][0]!=-1 and edges[i][1]!=-1){
                adj_list[edges[i][0]].push_back(edges[i][1]);
                adj_list[edges[i][1]].push_back(edges[i][0]);
                v1=edges[i][0];
                v2=edges[i][1];
                bool add_v1=true;
                bool add_v2=true;
                for(unsigned int i=0; i<nodes.size();i++){
                    if (nodes[i]==v1){
                        add_v1=false;
                    }
                    if (nodes[i]==v2){
                        add_v2=false;
                    }
                }
                if (add_v1==true)
                    nodes.push_back(v1);
                if (add_v2==true)
                    nodes.push_back(v2);
            }
        }

        for (unsigned int i=0;i<nodes.size();i++){
            neighbours[nodes[i]][0]=nodes[i];
            neighbours[nodes[i]][1]=adj_list[nodes[i]].size();
        }

        for (unsigned int i=0;i<nodes.size();i++){
            if (maximum<neighbours[nodes[i]][1] and neighbours[nodes[i]][1]!=-1){
                maximum=neighbours[nodes[i]][1];
                vertex=nodes[i];
             }//end if
        }//end for

        approxvc1_result.push_back(vertex);
		this->apvc1res.push_back(vertex);

        for (unsigned int i=0;i<no_edges;i++){
            if (edges[i][0]==vertex or edges[i][1]==vertex){
                edges[i][0]=-1;
                edges[i][1]=-1;

            }
        }
        all=false;
        for (unsigned int i=0;i<no_edges;i++){
            if(edges[i][0]!=-1){
                all=true;
                break;
            }
        }
    }
    while(all);

    sort(approxvc1_result.begin(),approxvc1_result.end());
	// std::cout<<"APPROX-VC-1: ";
    // for(unsigned int z=0;z<approxvc1_result.size();z++){
    //     std::cout<<approxvc1_result[z]<<" ";
    // }
    // cout<<endl;
}

// --------------Logic for approx-vc-2-----------------------//
int Param_Class_::approx_vc_2_recurrsivefunc(std::vector<int>approxvc2_list, int edgepair){
	//int smallno;
	std::vector<int> sortedlist;
	int allelem=0;

	//This loop will find the smallest no. in the input list. That no. will be the starting point
	for(int i=0;i<edgepair*2; i++){
		for (int j=0 ; j < (edgepair*2) ; j++){
			if(i==approxvc2_list[j] && approxvc2_list[j]!=-1){
				//smallno=approxvc2_list[j];
				sortedlist.push_back(approxvc2_list[j]);
				if(j%2==0)
					sortedlist.push_back(approxvc2_list[j+1]);
				else
				{
					sortedlist.push_back(approxvc2_list[j-1]);

				}
			}
		}
	}
	approxvc2_result.push_back(sortedlist[0]);
	approxvc2_result.push_back(sortedlist[1]);

	this->apvc2res.push_back(sortedlist[0]);
	this->apvc2res.push_back(sortedlist[1]);


	for(unsigned int i=0;i<approxvc2_list.size(); i++){
		if((approxvc2_list[i]==sortedlist[0] || approxvc2_list[i]==sortedlist[1]) && i%2==0)
		{
			approxvc2_list[i]=-1;
			approxvc2_list[i+1]=-1;
		}
		else if((approxvc2_list[i]==sortedlist[0] || approxvc2_list[i]==sortedlist[1]) && i%2!=0){
			approxvc2_list[i]=-1;
			approxvc2_list[i-1]=-1;

		}
	}
	sortedlist.clear();

	for(int i=0;i<edgepair*2; i++){
			if(approxvc2_list[i]==-1){
				allelem++;
			}
	}

	if(allelem!=edgepair*2){
		this->approx_vc_2_recurrsivefunc(approxvc2_list,edgepair);
	}
	return 0;
}

int Param_Class_::APPROX_VC_2(int edgepair){
	std::vector<int> approxvc2_list;

	for(int k=0 ; k<edgepair*2; k++){
		approxvc2_list.push_back(temp1dedges[k]);

	}

	this->approx_vc_2_recurrsivefunc(approxvc2_list,edgepair);
	// std::cout<<"APPROX-VC-2: ";
	// for(unsigned int z=0;z<approxvc2_result.size();z++){
	// 	std::cout<<approxvc2_result[z]<<" ";
	// }
	// std::cout<<"\n";

		// std::cout<<"Here  in apprx vc 2 function "<<!vertexcoverresult.empty()<< !vertexcoverfinalresult.empty()<<!approxvc1_result.empty()<<!approxvc2_result.empty()<<std::endl;
    return 0;
}

// --------------Logic for cnf-sat-----------------------//
std::vector<int> Param_Class_::get_vertexCover(vector<vector<int>> adj_list, int V) {
    //auto t_start=std::chrono::high_resolution_clock::now();
    if(adj_list.size()==0){
            std::cout<<endl;
			exit(0);
	    }


    for(int k=1;k<V;k++) {//parent for

        std::unique_ptr<Minisat::Solver> solver(new Minisat::Solver());
        vector<vector<Minisat::Lit>> expn(V);

        for(int i=0;i<V;i++){
            for(int j=0;j<k;j++) {
                Minisat::Lit l = Minisat::mkLit(solver->newVar());
                expn[i].push_back(l);
            }
        }

        //Clause 1 : At least one vertex is the ith vertex in the vertex cover
        for(int i=0;i<k;i++) {

            Minisat::vec<Minisat::Lit> clause;
            for(int n=0;n<V;n++) {
                clause.push(expn[n][i]);
            }
            solver->addClause(clause);
        }

        //Clause 2: No one vertex can appear twice in a vertex cover
        for(int m=0;m<V;m++){

            for(int p=0;p<k;p++){
                for(int q=0;q<k;q++){
                    if(p<q){
                        solver->addClause(~expn[m][p],~expn[m][q]);
                    }
                }
            }
        }

        //Clause 3: No more than one vertex appears in the mth position of the vertex cover
        for(int m=0;m<k;m++){
            for(int p=0;p<V;p++){
                for(int q=0; q<V; q++){
                    if(p<q){
                        solver->addClause(~expn[p][m],~expn[q][m]);
                    }
                }
            }
        }

        //Clause 4: Every edge is incident to at least one vertex in the vertex cover.
        for(unsigned int i = 0; i < adj_list.size(); i++) {
            Minisat::vec<Minisat::Lit> clause;
            for(int j = 0; j < k; j++) {
                clause.push(expn[adj_list[i][0]][j]);
                clause.push(expn[adj_list[i][1]][j]);
            }
            solver->addClause(clause);
        }

        //check if it is SAT and print output


        if(solver->solve()) {
			// std::cout<<"MINISAT: ";

            for(int n = 0; n < V; n++){
                for(int i = 0; i < k; i++){
                    if(solver->modelValue(expn[n][i]) == Minisat::l_True){
						vertexcoverresult.push_back(n);
                        // std::cout<<n<<" ";
						statusflag=200;

                    }
                }
            }
            // std::cout<<endl;
            solver.reset(new Minisat::Solver());

            return vertexcoverresult;
        }
        /*auto t_stop=std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t_stop-t_start;
		std::cout<<"Elapsed time: "<<elapsed.count()<<" s\n";
		if (elapsed.count()>=120){
            timeout=1;
            break;
		}*/
        solver.reset(new Minisat::Solver());
    }//end parent for

    return vertexcoverresult;

}
vector<vector<int>> adj_list_temp;

void Param_Class_::get_input(string line, int &V){

		pthread_mutex_lock( &mutex1 );

		// for(int i=0; i<edgepair*2;i++){
		// 	if(i%2==0)
		// 		adj_list.push_back({temp1dedges[i],temp1dedges[i+1]});
		// }


        this->vcres= get_vertexCover(adj_list_temp,Vertices);
		pthread_mutex_unlock( &mutex1 );

}

// --------------Function to print the final result-----------------------//

void Param_Class_::print(){
     if(!this->vcres.empty() && !this->apvc1res.empty() && !this->apvc2res.empty()){
			// std::cout<<"in print() vcres"<<!this->vcres.empty() <<std::endl;

        sort(this->vcres.begin(),this->vcres.end());
        sort(this->apvc1res.begin(),this->apvc1res.end());
        sort(this->apvc2res.begin(),this->apvc2res.end());

        std::cout<<"CNF-SAT-VC: ";
        for(unsigned int z=0;z<this->vcres.size();z++){
            std::cout<<this->vcres[z];
            if (z!=this->vcres.size()-1){
                std::cout<<",";
            }
        }
        cout<<endl;

        // std::cout<<"in print() apvc1res"<< !this->apvc1res.empty()<<std::endl;
        std::cout<<"APPROX-VC-1: ";
        for(unsigned int z=0;z<this->apvc1res.size();z++){
            std::cout<<this->apvc1res[z];
            if (z!=this->apvc1res.size()-1){
                std::cout<<",";
            }
        }
        cout<<endl;

        std::cout<<"APPROX-VC-2: ";
        for(unsigned int z=0;z<this->apvc2res.size();z++){
            std::cout<<this->apvc2res[z];
            if (z!=this->apvc2res.size()-1){
                std::cout<<",";
            }
        }
        std::cout<<"\n";
        result_flag=1;

        this->vcres.clear();
        this->apvc1res.clear();
        this->apvc2res.clear();

	}

	// else if(!this->vcres.empty() && !this->apvc1res.empty())
	// {
	// 		std::cout<<"in print() vcres"<<!this->vcres.empty() <<std::endl;

	// 			std::cout<<"in print() apvc1res"<< !this->apvc1res.empty()<<std::endl;

	// }
}

// --------------Threads-----------------------//
void* IO_thread(void *threaddata){
    //cout<<"in i/o thread"<<endl;
	//int starting_point;
	//int end_point;
	string Errorflag;
	int breakstatus=0;
    int start_pos=0, end_pos=0,n=0,i=0;
    int v_flag=0;
    int e_flag=0;
    vector <int> temp;


		// std::cout<<"Here  in io "<<!vertexcoverresult.empty()<<" "<< !vertexcoverfinalresult.empty()<<std::endl;
		while(!std::cin.eof() && breakstatus==0){

			std::getline(std::cin, inputval);
			if(std::cin.eof()){
                eof_all=1;
			}
			temp1dedges.clear();

			int spacecheck = inputval.find(' ');
			if (spacecheck > 0 )
			{
				if (inputval[0] == 'V') {
                    v_flag=1;

					// std::vector<std::string> splittedStrings = splitfunction(inputval, ' ');
					// stringstream NoOfVertices(splittedStrings[1]);
					// NoOfVertices >> Vertices;
					Vertices=stoi(inputval.substr(2,inputval.size()));


				}
				else if (inputval[0] == 'E')
				{
				    e_flag=1;
					// std::vector<std::string> splittedStrings = splitfunction(inputval, ' ');
					// TempEdges = splittedStrings[1];

					// for (int y = 0; y< TempEdges.length(); y++)
					// {
					// 	if (TempEdges[y] == '{')
					// 		TempEdges[y] = ' ';
					// 	else if (TempEdges[y] == '<')
                    //     {
                    //         edgepair++;
                    //         TempEdges[y] = ' ';
                    //     }
					// 	else if (TempEdges[y] == '>')
					// 		TempEdges[y] = ' ';
					// 	else if (TempEdges[y] == '}')
					// 		TempEdges[y] = ' ';
					// 	else if (TempEdges[y] == ',')
					// 		TempEdges[y] = ' ';
					// }

					// string str = TempEdges;
					// extractIntegersFromWords(str);
					unsigned int index=0;
					while (index<inputval.size()){
											// std::cout<<"Inside of while 1"<<std::endl;

						if(inputval.find("<",index)<inputval.size()){
							i=inputval.find("<",index);
							index=i;
							n=n+1;
						}//end if
						index=index+1;
					}//end while
					index=0;
					vector <int> edges;
					while(index<inputval.size() && n>0){
																	// std::cout<<"Inside of while 2"<<std::endl;

						int v1=-1, v2=-1;
						start_pos=inputval.find("<",end_pos);
						end_pos=inputval.find(">",start_pos);
						int len_sub=end_pos-start_pos-1;
						string sub=inputval.substr(start_pos+1,len_sub);
						//cout << "substring: "<< sub << '\n';
						int com= sub.find(",");
						v1 = stoi(sub.substr(0,com));
						v2 = stoi(sub.substr(com+1,sub.size()));
						if(v1 >= Vertices || v1<0){
								cerr << "Error: Invalid vertice" << '\n';
								exit(1);
						}
						else if(v2 >= Vertices || v2 < 0){
								cerr << "Error: Invalid vertice" << '\n';
								exit(1);
						}
						else if(v1 == v2){
						cerr << "Error: Invalid vertice" << '\n';
						exit(1);
						}
						adj_list_temp.push_back({v1,v2});
						temp1dedges.push_back(v1);
						temp1dedges.push_back(v2);


						// cout << "v1: " << v1 << "\t;v2: " << v2 << '\n';
						index=index+1;
						n=n-1;
					}//end while
					edgepair=adj_list_temp.size();
                // for (int i = 0; i <= edgepair*2; i++)
					// {
					// 	if (temp1dedges[i] >= Vertices)
					// 		errorflagstatus = "True";
					// }

					// if (errorflagstatus == "True") {
					// 	cout << "Error: Vertex does not exist OR Invalid Vertex" << endl;
					// 	TempEdges = " ";
					// 	adj1->clear();
					// 	for (int arrin = 0; arrin <= 50; arrin++) {
					// 		temp1dedges[arrin] = 0;
					// 	}
					// }


				}
				else{
					std::cout<<"Invalid Input"<<std::endl;
				}

			}
			else
			{
				continue;
			}

			if (v_flag==1 && e_flag==1){
                if (Vertices==0 || temp1dedges.empty()){
                    cout<<"CNF-SAT-VC: "<<endl<<"APPROX-VC-1: "<<endl<<"APPROX-VC-2: "<<endl;
                }
			}

			if(Vertices!=0 && !temp1dedges.empty()){
				// std::cout<<"break executed"<<statusflag<<std::endl;
				breakstatus=1;
				break;
			}

		}
		return NULL;
}

void* cnfsat_thread_(void *threaddata){
    //cout<<"in cnf thread"<<endl;
    struct _threadData_ *cnf_sat_ = (struct _threadData_ *)threaddata;

	int V;
	if(Vertices==0 && temp1dedges.empty()){
		// if(!vertexcoverresult.empty()){
		// std::cout<<edgepair*2<<std::endl;
		return NULL;
	}
	// int success;

    // success=pthread_getcpuclockid(pthread_self(), &threadClockId1);
    // if (success!=0){
    //     std::cerr<<"Error: pthread_getcpuclockid "<<success<<endl;
    // }
    // clock_gettime(threadClockId1, &time_s1);
    //cout<<"start "<<long(time_s1.tv_sec * 1000000 + time_s1.tv_nsec / 1000)<<endl;
    //printf("start  %.9f\n",time_s1.tv_sec+(double)(time_s1.tv_nsec)/1000000000);

	cnf_sat_->inputgraph.get_input(inputval, V);


	// success=pthread_getcpuclockid(pthread_self(), &threadClockId1);
    // if (success!=0){
    //     std::cerr<<"Error: pthread_getcpuclockid "<<success<<endl;
    // }
    // clock_gettime(threadClockId1, &time_e1);

    // //cout<<"end"<<long(time_s2.tv_sec * 1000000 + time_s2.tv_nsec / 1000)<<endl;
    // //printf("end  %.9f\n",time_e1.tv_sec+(double)(time_e1.tv_nsec)/1000000000);
    // printf("cnf sat time: %.12f us\n", ((time_e1.tv_sec+(double)(time_e1.tv_nsec)/1000000000)-(time_s1.tv_sec+(double)(time_s1.tv_nsec)/1000000000))*1000000);

    cnf_sat_->inputgraph.print();
	return NULL;
}

void* approx_vc_1_thread_(void *threaddata){
			// std::cout<<"Here  in apprx vc 1 "<<!vertexcoverresult.empty()<<" "<< !vertexcoverfinalresult.empty()<<std::endl;
    //cout<<"in approx 1 thread"<<endl;
    struct _threadData_ *approxvc1_ = (struct _threadData_ *)threaddata;

	if(Vertices==0 && temp1dedges.empty()){
		// if(!vertexcoverresult.empty()){
		// std::cout<<edgepair*2<<std::endl;
		return NULL;
	}
	// int success;
	// success=pthread_getcpuclockid(pthread_self(), &threadClockId2);
    // if (success!=0){
    //     std::cerr<<"Error: pthread_getcpuclockid "<<success<<endl;
    // }
    // clock_gettime(threadClockId2, &time_s2);


	approxvc1_->inputgraph.APPROX_VC_1(Vertices,edgepair);


	// success=pthread_getcpuclockid(pthread_self(), &threadClockId2);
    // if (success!=0){
    //     std::cerr<<"Error: pthread_getcpuclockid "<<success<<endl;
    // }
    // clock_gettime(threadClockId2, &time_e2);

    // printf("vc1 time: %.12f us\n", ((time_e2.tv_sec+(double)(time_e2.tv_nsec)/1000000000)-(time_s2.tv_sec+(double)(time_s2.tv_nsec)/1000000000))*1000000);

    approxvc1_->inputgraph.print();
	return NULL;
}

void* approx_vc_2_thread_(void *threaddata){
    //cout<<"in approx 2 thread"<<endl;
	struct _threadData_ *approxvc2_ = (struct _threadData_ *)threaddata;

	if(Vertices==0 && temp1dedges.empty()){
		return NULL;
	}
	// int success;

	// success=pthread_getcpuclockid(pthread_self(), &threadClockId3);
    // if (success!=0){
    //     std::cerr<<"Error: pthread_getcpuclockid "<<success<<endl;
    // }
    // clock_gettime(threadClockId3, &time_s3);

	approxvc2_->inputgraph.APPROX_VC_2(edgepair);


	// success=pthread_getcpuclockid(pthread_self(), &threadClockId3);
    // if (success!=0){
    //     std::cerr<<"Error: pthread_getcpuclockid "<<success<<endl;
    // }
    // clock_gettime(threadClockId3, &time_e3);

    // printf("vc 2 time: %.12f us\n", ((time_e3.tv_sec+(double)(time_e3.tv_nsec)/1000000000)-(time_s3.tv_sec+(double)(time_s3.tv_nsec)/1000000000))*1000000);

    approxvc2_->inputgraph.print();
	return NULL;
}

/*void* printresult(void *threaddata){
	// std::cout<<"In printresultfunc" <<" "<<!vertexcoverfinalresult.empty() << !approxvc1_result.empty()<<!approxvc2_result.empty()<<std::endl;

	if(!vertexcoverfinalresult.empty() && !approxvc1_result.empty() && !approxvc2_result.empty())
		{
 			// std::cout<<"MINISAT: ";
			// for (int i=0; i<vertexcoverfinalresult.size();i++){
			// 	std::cout<<vertexcoverfinalresult[i]<<" ";
			// }

			std::cout<<"APPROX-VC-1: n";
			for(unsigned int z=0;z<approxvc1_result.size();z++){
				std::cout<<approxvc1_result[z]<<" ";
			}
			std::cout<<"\n";

			std::cout<<"APPROX-VC-2: ";
			for(int i=0; i<approxvc2_result.size(); i++){
				std::cout<<approxvc2_result[i]<<" ";
			}
			std::cout<<"\n";

			vertexcoverfinalresult.clear();
			approxvc2_result.clear();
			approxvc2_result.clear();
		}

	return NULL;
}*/

int main(int argc, char** argv)
{


    struct _threadData_ inputgraph;
	pthread_t io_thread;
    pthread_t cnfsat_thread;
    pthread_t approx_vc_1_thread;
	pthread_t approx_vc_2_thread;

    auto t_start=std::chrono::high_resolution_clock::now();
    int l_count=0;


	while(true){
        l_count=l_count+1;
        //cout<<"in while"<<approxvc1_result.size()<<" "<<approxvc2_result.size()<<endl;
		// srand(time(nullptr)); // srand() in main functoin to ensure the unique random value of each time
		//  int ret;
    	// for(size_t i = 0; i < 5;++i) {
		// 	ret = pthread_create(&threads[i], nullptr, algorithm[i], &statusflag);
		// 	if (ret) _exit(1); /*error happens here*/
    	// }

    	if (eof_all == 1){
            break;
    	}


		pthread_create(&io_thread, NULL, IO_thread, &inputgraph);
		pthread_create(&cnfsat_thread, NULL, cnfsat_thread_, &inputgraph);
		pthread_create(&approx_vc_1_thread, NULL, approx_vc_1_thread_, &inputgraph);
		pthread_create(&approx_vc_2_thread, NULL, approx_vc_2_thread_, &inputgraph);

        if (l_count!=1){
            while(timeout==0 and result_flag==0){
                auto t_stop=std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = t_stop-t_start;
                //std::cout<<"Elapsed time: "<<elapsed.count()<<" s\n";
                if (elapsed.count()>=120){
                    timeout=1;
                }
            }
        }
        if (timeout==1){
            break;
        }

		// for(size_t i = 0; i < 5;++i) {
        // pthread_join(threads[i],NULL);
		// }
		pthread_join(io_thread, NULL);
		pthread_join(cnfsat_thread, NULL);
		pthread_join(approx_vc_1_thread, NULL);
		pthread_join(approx_vc_2_thread, NULL);




	}
	//cout<<"exited"<<endl;

	sort(approxvc1_result.begin(),approxvc1_result.end());
	sort(approxvc2_result.begin(),approxvc2_result.end());

	if (timeout==1){
        std::cout<<"CNF-SAT-VC: timeout";
        cout<<endl;

        // std::cout<<"in print() apvc1res"<< !this->apvc1res.empty()<<std::endl;
        std::cout<<"APPROX-VC-1: ";
        for(unsigned int z=0;z<approxvc1_result.size();z++){
            std::cout<<approxvc1_result[z];
            if (z!=approxvc1_result.size()-1){
                std::cout<<",";
            }
        }
        cout<<endl;

        std::cout<<"APPROX-VC-2: ";
        for(unsigned int z=0;z<approxvc2_result.size();z++){
            std::cout<<approxvc2_result[z];
            if (z!=approxvc2_result.size()-1){
                std::cout<<",";
            }
        }
        std::cout<<"\n";

	}

	return 0;


}
