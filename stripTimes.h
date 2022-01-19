#ifndef stripTimes_h
#define stripTimes_h


using namespace std;

struct stripTimes
{
        void addLRframe(int i , float time) {LRframe.push_back(i); LRtime.push_back(time) ; LRframeOK.push_back(true); LRframeBkg.push_back(false);}
        void addHRframe(int i, float time) {HRframe.push_back(i); HRtime.push_back(time) ; HRframeOK.push_back(true); HRframeBkg.push_back(false);}
        std::vector<int> LRframe;
        std::vector<float> LRtime;
        std::vector<bool> LRframeOK;
        std::vector<int> HRframe;
        std::vector<float> HRtime;
        std::vector<bool> HRframeOK;
        std::vector<bool> HRframeBkg;
        std::vector<bool> LRframeBkg;

        std::vector<std::pair<int,int> > deltas;
        void buildDeltas()
        {
                deltas.clear();
                for (int iLR=0; iLR<LRframe.size();++iLR)
                {
                        if ( not (LRframeOK[iLR])) continue;
                        for (int iHR=0; iHR<HRframe.size(); ++iHR)
                        {
                                if (HRframeOK[iHR]) deltas.push_back(std::make_pair(LRframe[iLR],HRframe[iHR]));
                        }
                }
        }
};

#endif
