#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <utility>

struct EtaBin {
  std::pair<float, float> bin;
  std::string name;
  std::string title;
};

class EtaBinning {
  public:
    EtaBinning() {
      fillEtaBins();
    }

    int getBin(float eta) const {
      eta = fabs(eta);
      std::vector<EtaBin>::const_iterator it = mEtaBins.begin();
      for (; it != mEtaBins.end(); ++it) {
        EtaBin bin = *it;
        if (eta >= bin.bin.first && eta < bin.bin.second) {
          return it - mEtaBins.begin();
        }
      }

      return -1;
    }

    std::string getBinName(int bin) const {
      return mEtaBins[bin].name;
    }

    std::string getBinTitle(int bin) const {
      return mEtaBins[bin].title;
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mEtaBins[bin].bin;
    }


    size_t size() const {
      return mEtaBins.size();
    }

  private:
    std::vector<EtaBin> mEtaBins;

    void fillEtaBins() {

      EtaBin bin;

      bin.bin = std::make_pair(0., 0.5);
      bin.name = "eta_00_05";
      bin.title = "|#eta| < 0.5";  
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(0.5, 1.0);
      bin.name = "eta_05_10";
      bin.title = "0.5 #leq |#eta| < 1.0";  
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(1.0, 1.5);
      bin.name = "eta_10_15";
      bin.title = "1.3 #leq |#eta| < 2.0";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.5, 2.0);
      bin.name = "eta_15_20";
      bin.title = "1.5 #leq |#eta| < 2.0";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(2.0, 2.5);
      bin.name = "eta_20_25";
      bin.title = "2.0 #leq |#eta| < 2.5";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(2.5, 3.0);
      bin.name = "eta_25_30";
      bin.title = "2.5 #leq |#eta| < 3.0";
      mEtaBins.push_back(bin);


    }
};
