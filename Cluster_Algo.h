#ifndef Cluster_Algo_h
#define Cluster_Algo_h


using namespace std;

struct cluster{
        void addClusterPair(std::pair<vector<int>,vector<double>> LR_v, std::pair<vector<int>,vector<double>> HR_v, double shift_deltaT_){strip_LR=LR_v.first; time_LR=LR_v.second; strip_HR=HR_v.first; time_HR=HR_v.second; shift_deltaT=shift_deltaT_;}
        std::vector<int> strip_LR;
        std::vector<double> time_LR;
        std::vector<int> strip_HR;
        std::vector<double> time_HR;
	double shift_deltaT;

        double time(){
                double timeLR = 0.;
                for(int i=0;i<time_LR.size();i++) {timeLR=timeLR+time_LR.at(i);}
                timeLR=timeLR/time_LR.size();
                double timeHR = 0.;
                for(int i=0;i<time_HR.size();i++) {timeHR=timeHR+time_HR.at(i);}
                timeHR=timeHR/time_HR.size();

                double time = timeLR-timeHR+shift_deltaT;
		time = time/2.;
                return time;
        }

        double sum_time(){
                double timeLR = 0.;
                for(int i=0;i<time_LR.size();i++) {timeLR=timeLR+time_LR.at(i);}
                timeLR=timeLR/time_LR.size();
                double timeHR = 0.;
                for(int i=0;i<time_HR.size();i++) {timeHR=timeHR+time_HR.at(i);}
                timeHR=timeHR/time_HR.size();

                double time = timeLR+timeHR;
                return time;
        }



        double strip(){
                double strip_number = ((double)std::accumulate(strip_LR.begin(), strip_LR.end(), 0)+(double)std::accumulate(strip_HR.begin(), strip_HR.end(), 0))/(strip_LR.size()+strip_HR.size());
                return strip_number;
        }

        int size(){
                std::vector<int> intersec_vec_temp;
                intersec_vec_temp.clear();
                set_intersection(strip_LR.begin(),strip_LR.end(),strip_HR.begin(),strip_HR.end(),back_inserter(intersec_vec_temp));
                int size = strip_HR.size()+strip_LR.size()-intersec_vec_temp.size();
                return size;
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



vector<std::pair<vector<int>,vector<double>>> cluster_one_side(vector<int> strip, vector<double> time, double limit){

        vector<std::pair<vector<int>,vector<double>>> cluster_one_side;
        cluster_one_side.clear();
	
	vector<int> indice;
	indice.clear();
//	while(strip.size()>0){
 	       for(int j=0; j<strip.size(); j++){
        	        vector<int> v_strip_temp;
                	vector<double> v_time_temp;
                	v_strip_temp.clear();
                	v_time_temp.clear();
                	v_strip_temp.push_back(strip.at(j));
                	v_time_temp.push_back(time.at(j));
			indice.push_back(j);
                	if((j-strip.size())==-1){ // if we are in the last strip
                        	cluster_one_side.push_back(std::make_pair(v_strip_temp,v_time_temp));
                        	break;
                	} else {
                        	for(int i=j+1; i<strip.size(); i++){
                                	if(abs(time.at(i-1)-time.at(i))<limit && ((strip.at(i-1)-strip.at(i))==-1 || (strip.at(i-1)-strip.at(i))==-2)){
                                        	v_strip_temp.push_back(strip.at(i));
                                        	v_time_temp.push_back(time.at(i));
						indice.push_back(i);
					} else if((i+1)<strip.size() && abs(time.at(i-1)-time.at(i+1))<limit && (strip.at(i-1)-strip.at(i+1))==-2){
                                        	v_strip_temp.push_back(strip.at(i+1));
                                        	v_time_temp.push_back(time.at(i+1));
						indice.push_back(i+1);
						i++;
                                	} else {
                                        	cluster_one_side.push_back(std::make_pair(v_strip_temp,v_time_temp));
                                        	j = i-1;
                                        	continue;
                                	}
                        	}
                	}
        	}
//                for(int j=(indice.size()-1);j>=0;j--){
//			strip.erase(strip.begin()+indice.at(j));
//			time.erase(time.begin()+indice.at(j));
//		}
//	}
        return cluster_one_side;
}



std::vector<cluster> final_clustering(vector<std::pair<vector<int>,vector<double>>> LR_c,vector<std::pair<vector<int>,vector<double>>> HR_c, double shift_deltaT){

 	bool print=false;
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




#endif
