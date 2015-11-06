#ifndef FWLite_TFileService_h
#define FWLite_TFileService_h
#include "TFileDirectory.h"

namespace fwlite {
  
  class TFileService : public TFileDirectory {
  public:
    TFileService(const std::string& fileName);
    
    TFileService(TFile * aFile);
    
    ~TFileService();
    
    TFile & file() const { return * file_; }
    
  private:
    TFile * file_;
    std::string fileName_;
    
  };
  
}
#endif

