#pragma once

#include <cmath>
#include <vector>
#include <utility>

class PtBinning {
  public:
    PtBinning() {
      fillPtBins();
    }

    int getPtBin(float pt) {
      std::vector<std::pair<float, float> >::const_iterator it = mPtBins.begin();
      for (; it != mPtBins.end(); ++it) {
        std::pair<float, float> bin = *it;
        if (pt >= bin.first && pt < bin.second) {
          return it - mPtBins.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mPtBins.size();
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mPtBins[bin];
    }

    std::vector<std::pair<float, float> > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<std::pair<float, float> >(mPtBins.begin(), mPtBins.begin() + n);
    }

    std::vector<std::pair<float, float> > getBinning(unsigned int from, unsigned int to) const {
      if (to > size()) {
        to = size();
      }

      return std::vector<std::pair<float, float> >(mPtBins.begin() + from, mPtBins.begin() + to);
    }

    int get_PtStep(){   
      return step_pt;
    }


  private:
    std::vector<std::pair<float, float> > mPtBins;
    int step_pt = 200;
    int pT_max= 5000;

    void fillPtBins() {
      
      int n_bin = pT_max/ step_pt;

      for(int jj=0; jj<n_bin ; jj++){	 
	int pt_min = ( (jj-1)*step_pt)+step_pt;
	int pt_max = (jj*step_pt)+step_pt;
	
      mPtBins.push_back(std::make_pair(pt_min, pt_max));

	}
    } 
};

      //      mPtBins.push_back(std::make_pair(60., 85.));
      //mPtBins.push_back(std::make_pair(85., 100.));
      // mPtBins.push_back(std::make_pair(100., 130.));
      // mPtBins.push_back(std::make_pair(130., 175.));
      // mPtBins.push_back(std::make_pair(175., 250.));
      // mPtBins.push_back(std::make_pair(250., 300.));
      // mPtBins.push_back(std::make_pair(300., 400.));
      // mPtBins.push_back(std::make_pair(400., 1100.));
