#ifndef Cluster_Algo_h
#define Cluster_Algo_h
#include "CalibAlig.h"


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

        vector<double> time_v(){
                return time_LR;
        }
};


struct cluster{
	void addClusterPair(ClusterOneSide LR_v, ClusterOneSide HR_v, double shift_deltaT_, int sn__){strip_LR=LR_v.strip(); time_LR=LR_v.time(); strip_HR=HR_v.strip(); time_HR=HR_v.time(); size_LR=LR_v.size(); size_HR=HR_v.size(); strip_LR_v=LR_v.strips_v(); strip_HR_v=HR_v.strips_v(); shift_deltaT=shift_deltaT_; sn_=sn__; time_LR_v=LR_v.time_v(); time_HR_v=HR_v.time_v();} 

	int sn_;
	double strip_LR;
	double time_LR;
	double strip_HR;
	double time_HR;
	double shift_deltaT;
        int size_LR;
        int size_HR;
	vector<int> strip_LR_v;
        vector<int> strip_HR_v;
        vector<double> time_LR_v;
        vector<double> time_HR_v;

	double time(){
		double time = time_LR-time_HR+shift_deltaT;
		time = time/2.;
		return time;
	}

        double time_XY(){
                double time = time_HR-time_LR;
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

	double X(){
		double x = 0.;
		int n = 0;
		for(int i=0;i<strip_HR_v.size();i++){
			for(int j=0;j<strip_LR_v.size();j++){
				if(strip_LR_v.at(j)==strip_HR_v.at(i)){
                                        vector<double> X_Y = convert_DeltaTAndStrip_To_XY(sn_,strip_HR_v.at(i),(time_HR_v.at(i)-time_LR_v.at(j)-dT_align_factor[strip_LR_v.at(j)]));
					x = x + X_Y.at(0);
					n = n + 1;
				}
			}
		}

		if(n!=0) x = x/n;
		return x;
	}

	double Y(){
                double y = 0.;
                int n = 0;
                for(int i=0;i<strip_HR_v.size();i++){
                        for(int j=0;j<strip_LR_v.size();j++){
                                if(strip_LR_v.at(j)==strip_HR_v.at(i)){
                                        vector<double> X_Y = convert_DeltaTAndStrip_To_XY(sn_,strip_HR_v.at(i),(time_HR_v.at(i)-time_LR_v.at(j)-dT_align_factor[strip_LR_v.at(j)]));
                                        y = y + X_Y.at(1);
                                        n = n + 1;
                                }
                        }
                }

                if(n!=0) y = y/n;
		return y;
	}

	double HR_LR_strips_ratio(){
		double ratio = (double)strip_HR_v.size()/(double)strip_LR_v.size();
		return ratio;
	}

};

vector<ClusterOneSide> cluster_one_side(vector<int> strip, vector<double> time, double limit){

	vector<std::pair<vector<int>,vector<double>>> cluster_one_side;
	//	vector<ClusterOneSide> cluster_one_side;
	cluster_one_side.clear();

	vector<int> indice;
	indice.clear();
	for(int j=0; j<time.size(); j++){
		indice.push_back(j);
	}


	while(indice.size()>0){
		double max_time = 9999999999999999999999999999.;
		int max_time_ind = -999;
		for(int j=0; j<indice.size(); j++){
			if(time.at(indice.at(j))<max_time) {max_time = time.at(indice.at(j)); max_time_ind = j;}
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
		int time_ref =-99;
		for(int j=max_time_ind; j<indice.size(); j++){
			if(j==max_time_ind) {i=max_time_ind; time_ref=max_time_ind; continue;}
//                        if(abs(time.at(indice.at(i))-time.at(indice.at(j)))<limit){
			if(abs(time.at(indice.at(time_ref))-time.at(indice.at(j)))<limit){ 
				if(abs(strip.at(indice.at(i))-strip.at(indice.at(j)))<=1+max_jumps){
					v_strip_temp.push_back(strip.at(indice.at(j)));
					v_time_temp.push_back(time.at(indice.at(j)));
					indice_to_erase.push_back(j);
					i=j;
				} else break;
			}
		}

		i=-99;
		time_ref=-99;
		for(int j=max_time_ind; j>-1; j--){
			if(j==max_time_ind) {i=max_time_ind; time_ref=max_time_ind; continue;}
                        if(abs(time.at(indice.at(time_ref))-time.at(indice.at(j)))<limit){
//			if(abs(time.at(indice.at(i))-time.at(indice.at(j)))<limit){
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

	return cluster_one_side_clust;
}


std::vector<cluster> final_clustering(vector<ClusterOneSide> LR_c, vector<ClusterOneSide> HR_c,double shift_deltaT, int sn_){

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
		if(abs(HR_c.at(HR_ind_smallest).strip() - LR_c.at(LR_ind_smallest).strip())<=1.){
			cluster temp0;
			temp0.addClusterPair(LR_c.at(LR_ind_smallest),HR_c.at(HR_ind_smallest),shift_deltaT,sn_);
			Cluster_v.push_back(temp0);
			LR_c.erase(LR_c.begin()+LR_ind_smallest);
			HR_c.erase(HR_c.begin()+HR_ind_smallest);
		} else {no_match++;}
	}

	return Cluster_v;
}

#endif
