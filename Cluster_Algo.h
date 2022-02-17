#ifndef Cluster_Algo_h
#define Cluster_Algo_h


using namespace std;

struct ClusterOneSide{
	void addClusterPair(std::pair<vector<int>,vector<double>> LR_v){strip_LR=LR_v.first; time_LR=LR_v.second;}
	std::vector<int> strip_LR;
	std::vector<double> time_LR;

	double time(){
		double timeLR = 0.;
		for(int i=0;i<time_LR.size();i++) {timeLR=timeLR+time_LR.at(i);}
		timeLR=timeLR/time_LR.size();
		return timeLR;
	}

	double strip(){
		double stripLR = 0.;
		for(int i=0;i<strip_LR.size();i++) {stripLR=stripLR+strip_LR.at(i);}
		stripLR=stripLR/strip_LR.size();
		return stripLR;
	}

	int size(){
		return time_LR.size();
	}

	vector<int> strips_v(){
		return strip_LR;
	}

};


struct cluster{
	void addClusterPair(ClusterOneSide LR_v, ClusterOneSide HR_v, double shift_deltaT_){strip_LR=LR_v.strip(); time_LR=LR_v.time(); strip_HR=HR_v.strip(); time_HR=HR_v.time(); size_LR=LR_v.size(); size_HR=HR_v.size(); strip_LR_v=LR_v.strips_v(); strip_HR_v=HR_v.strips_v(); shift_deltaT=shift_deltaT_;}
	double strip_LR;
	double time_LR;
	double strip_HR;
	double time_HR;
	double shift_deltaT;
        int size_LR;
        int size_HR;
	vector<int> strip_LR_v;
        vector<int> strip_HR_v;

	double time(){
		double time = time_LR-time_HR+shift_deltaT;
		time = time/2.;
		return time;
	}

	double sum_time(){
		double time = time_LR+time_HR;
		return time;
	}



	double strip(){
		double strip = strip_LR+strip_HR;
		strip = strip/2;
		return strip;
	}

	int size(){
		int size = 0;
		for(int i=0;i<strip_LR_v.size();i++){
        	        for(int j=0;j<strip_HR_v.size();j++){
				if(strip_LR_v.at(i)==strip_HR_v.at(j)) size++;
			}
		}
		size = -size+strip_LR_v.size()+strip_HR_v.size();
		return size;
	}

	double HR_LR_strips_ratio(){
		double ratio = (double)strip_HR_v.size()/(double)strip_LR_v.size();
		return ratio;
	}

};


struct cluster_new_try{
	void addClusterPair(std::vector<int> strip_, std::vector<double> time_){time_v=time_;strip_v=strip_;}
	std::vector<int> strip_v;
	std::vector<double> time_v;
	double time(){
		double dtime = 0.;
		for(int i=0;i<time_v.size();i++) {dtime=dtime+time_v.at(i);}
		dtime=dtime/time_v.size();
		return dtime;
	}
	double strip(){
		double dstrip = 0.;
		for(int i=0;i<strip_v.size();i++) {dstrip=dstrip+strip_v.at(i);}
		dstrip=dstrip/strip_v.size();
		return dstrip;
	}
	int size(){
		return strip_v.size();
	}
};


vector<cluster_new_try> cluster_new_try_func(vector<double> strip, vector<double> time){
	vector<cluster_new_try> CLUSTER;
	CLUSTER.clear();

	double limit = 1.;
	vector<double> strip_temp = strip;
	vector<double> time_temp = time;

	while(strip_temp.size()>0){
		vector<int> indices;
		indices.clear();

		vector<int> v_strip_temp;
		vector<double> v_time_temp;
		v_strip_temp.clear();
		v_time_temp.clear();


		v_strip_temp.push_back(strip.at(0));
		v_time_temp.push_back(time.at(0));
		indices.push_back(0);

		int cont = 0;
		int j = 0;
		for(int i = 1; i<strip_temp.size();i++){
			if(cont==2) break;
			if( abs(time_temp.at(i)-time.at(j))<limit && ((strip_temp.at(i)-strip_temp.at(j))==1)){
				v_strip_temp.push_back(strip.at(i));
				v_time_temp.push_back(time.at(i));
				indices.push_back(i);
				j = i;
				cont=0;
			} else  cont++;
		}

		cluster_new_try new_cluster;
		new_cluster.addClusterPair(v_strip_temp, v_time_temp);
		CLUSTER.push_back(new_cluster);
		for(int j=(indices.size()-1);j>=0;j--){
			strip_temp.erase(strip_temp.begin()+indices.at(j));
			time_temp.erase(time_temp.begin()+indices.at(j));
		}

	}
	return CLUSTER;
}

//vector<std::pair<vector<int>,vector<double>>> cluster_one_side(vector<int> strip, vector<double> time, double limit){
vector<ClusterOneSide> cluster_one_side(vector<int> strip, vector<double> time, double limit){

	vector<std::pair<vector<int>,vector<double>>> cluster_one_side;
	//	vector<ClusterOneSide> cluster_one_side;
	cluster_one_side.clear();

	vector<int> indice;
	indice.clear();
	//cout << "****************************************" << endl;
	for(int j=0; j<time.size(); j++){
	//	cout << "strip.at(j): " << strip.at(j) << endl;
		indice.push_back(j);
	}


	while(indice.size()>0){
		double max_time = -9999999999999999999999999999.;
		int max_time_ind = -999;
		for(int j=0; j<indice.size(); j++){
			if(time.at(indice.at(j))>max_time) {max_time = time.at(indice.at(j)); max_time_ind = j;}
		}

		vector<int> v_strip_temp;
		vector<double> v_time_temp;
		v_strip_temp.clear();
		v_time_temp.clear();
		v_strip_temp.push_back(strip.at(indice.at(max_time_ind)));
		v_time_temp.push_back(time.at(indice.at(max_time_ind)));



		vector<int> indice_to_erase;
		indice_to_erase.clear();
		indice_to_erase.push_back(max_time_ind);

		int max_jumps = 1;

		int i = -99;
		for(int j=max_time_ind; j<indice.size(); j++){
			if(j==max_time_ind) {i=max_time_ind; continue;}
			if(abs(time.at(indice.at(i))-time.at(indice.at(j)))<limit){ 
				if(abs(strip.at(indice.at(i))-strip.at(indice.at(j)))<=1+max_jumps){
					v_strip_temp.push_back(strip.at(indice.at(j)));
					v_time_temp.push_back(time.at(indice.at(j)));
					indice_to_erase.push_back(j);
					i=j;
				} else break;
			}
		}

		i=-99;
		for(int j=max_time_ind; j>-1; j--){
			if(j==max_time_ind) {i=max_time_ind; continue;}
			if(abs(time.at(indice.at(i))-time.at(indice.at(j)))<limit){
				if(abs(strip.at(indice.at(i))-strip.at(indice.at(j)))<=1+max_jumps){
					v_strip_temp.push_back(strip.at(indice.at(j)));
					v_time_temp.push_back(time.at(indice.at(j)));
					indice_to_erase.push_back(j);
					i=j;
				} else break;
			}
		}

		sort(indice_to_erase.begin(), indice_to_erase.end());

		for(int j=indice_to_erase.size()-1; j>-1; j--){
			indice.erase(indice.begin()+indice_to_erase.at(j));
		}

		cluster_one_side.push_back(std::make_pair(v_strip_temp,v_time_temp));


	}

	vector<ClusterOneSide> cluster_one_side_clust;
	cluster_one_side_clust.clear();

	for(int j=0; j<cluster_one_side.size(); j++){
		ClusterOneSide new_cluster;
		new_cluster.addClusterPair(cluster_one_side.at(j));
		cluster_one_side_clust.push_back(new_cluster);
	}

/*        for(int j=0; j<cluster_one_side_clust.size(); j++){
	cout << "================================================" << endl;
		for(int i=0;i<cluster_one_side_clust.at(j).strips_v().size();i++){
			cout << " cluster_one_side_clust.at(j).strips_v().at(i): " << cluster_one_side_clust.at(j).strips_v().at(i) << endl;
		}
	}
	cout << endl;
*/
	return cluster_one_side_clust;
}



//std::vector<cluster> final_clustering(vector<std::pair<vector<int>,vector<double>>> LR_c,vector<std::pair<vector<int>,vector<double>>> HR_c, double shift_deltaT){

std::vector<cluster> final_clustering(vector<ClusterOneSide> LR_c, vector<ClusterOneSide> HR_c,double shift_deltaT){

	std::vector<cluster> Cluster_v;
	Cluster_v.clear();

	int no_match = 0;
	while(LR_c.size()>0 && HR_c.size()>0 && no_match==0){
		double strip_diff = 999.;
		int LR_ind_smallest = -1;
		int HR_ind_smallest = -1;
		for(int ij=0;ij<HR_c.size();ij++){
			for(int ji=0;ji<LR_c.size();ji++){
				if(abs(HR_c.at(ij).strip() - LR_c.at(ji).strip())<strip_diff) {
					strip_diff = abs(HR_c.at(ij).strip() - LR_c.at(ji).strip());
					LR_ind_smallest = ji;
					HR_ind_smallest = ij;
				}
			}
		}
		if(abs(HR_c.at(HR_ind_smallest).strip() - LR_c.at(LR_ind_smallest).strip())<0.9){
			cluster temp0;
			temp0.addClusterPair(LR_c.at(LR_ind_smallest),HR_c.at(HR_ind_smallest),shift_deltaT);
			//if(temp0.time()>0.){
				Cluster_v.push_back(temp0);
				LR_c.erase(LR_c.begin()+LR_ind_smallest);
				HR_c.erase(HR_c.begin()+HR_ind_smallest);
			//}
		} else {no_match++;}
	}

	return Cluster_v;
}
/*	bool print=false;
	std::vector<cluster> Cluster_v;
	Cluster_v.clear();

	vector<std::pair<vector<int>,vector<double>>> LR_c_clone = LR_c;
	vector<std::pair<vector<int>,vector<double>>> HR_c_clone = HR_c;

// combine pairs

while(LR_c.size()>0 && HR_c.size()>0){
int LR_i = -999;
int HR_i = -999;
int n_intersec = 0;
for(int i=0; i<LR_c.size(); i++){
for(int j=0; j<HR_c.size(); j++){
vector<int> intersec_vec_tem;
intersec_vec_tem.clear();
set_intersection(LR_c.at(i).first.begin(),LR_c.at(i).first.end(),HR_c.at(j).first.begin(),HR_c.at(j).first.end(),back_inserter(intersec_vec_tem));
cluster temp0;
temp0.addClusterPair(LR_c.at(i),HR_c.at(j),shift_deltaT);
if(intersec_vec_tem.size()>n_intersec && temp0.time()>0.){
LR_i = i;
HR_i = j;
n_intersec = intersec_vec_tem.size();
}
}
}
if(LR_i!=-999 && HR_i!=-999){
cluster temp;
temp.addClusterPair(LR_c.at(LR_i),HR_c.at(HR_i),shift_deltaT);// save the cluster pair
//              if((temp.time()+shift_deltaT)<0.){
Cluster_v.push_back(temp);
LR_c.erase(LR_c.begin()+LR_i);
HR_c.erase(HR_c.begin()+HR_i);
//              }
} else break;
}


if(Cluster_v.size()==0){
while(LR_c.size()>0 && HR_c.size()>0){
int LR_i = -999;
int HR_i = -999;
for(int i=0; i<LR_c.size(); i++){
for(int j=0; j<HR_c.size(); j++){

cluster temp0;
temp0.addClusterPair(LR_c.at(i),HR_c.at(j),shift_deltaT);
if((abs(LR_c.at(i).first.at(0)-HR_c.at(j).first.at(HR_c.at(j).first.size()-1))==1||abs(HR_c.at(j).first.at(0)-LR_c.at(i).first.at(LR_c.at(i).first.size()-1))==1||abs(LR_c.at(i).first.at(0)-HR_c.at(j).first.at(0))==1||abs(LR_c.at(i).first.at(LR_c.at(i).first.size()-1)-HR_c.at(j).first.at(HR_c.at(j).first.size()-1))==1) && (temp0.time())>0.){
LR_i = i;
HR_i = j;
}
}
}
if(LR_i!=-999 && HR_i!=-999){
cluster temp;
temp.addClusterPair(LR_c.at(LR_i),HR_c.at(HR_i),shift_deltaT);// save the cluster pair
//                     if((temp.time()+shift_deltaT)<0.){
Cluster_v.push_back(temp);
LR_c.erase(LR_c.begin()+LR_i);
HR_c.erase(HR_c.begin()+HR_i);
//                      }
} else break;
}
}


if(print){
if(Cluster_v.size()==0){
cout << "----------->>>  Cluster_v.size(): " << Cluster_v.size() << "; LR_c_clone.size(): " << LR_c_clone.size() << "; HR_c_clone.size(): " << HR_c_clone.size() << endl;
for(int ii=0; ii<LR_c_clone.size(); ii++){
cout << "-->> LR_c_clone.at(i).first.size(): " << LR_c_clone.at(ii).first.size() << endl;
for(int jj=0; jj<LR_c_clone.at(ii).first.size(); jj++){
	cout << "LR_c_clone.at(ii).first.at(jj): " << LR_c_clone.at(ii).first.at(jj) << endl;
}
}
for(int ii=0; ii<HR_c_clone.size(); ii++){
	cout << "-->> HR_c_clone.at(i).first.size(): " << HR_c_clone.at(ii).first.size() << endl;
	for(int jj=0; jj<HR_c_clone.at(ii).first.size(); jj++){
		cout << "HR_c_clone.at(ii).first.at(jj): " << HR_c_clone.at(ii).first.at(jj) << endl;
	}
}
}
}

return Cluster_v;
}
*/



#endif
