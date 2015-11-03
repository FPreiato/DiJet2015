#ifndef CommonTools_Utils_TH1AddDirectorySentry_h
#define CommonTools_Utils_TH1AddDirectorySentry_h

class TH1AddDirectorySentry
{
  
    public:
  TH1AddDirectorySentry();
  ~TH1AddDirectorySentry();
  
  
 private:
  TH1AddDirectorySentry(const TH1AddDirectorySentry&);
  TH1AddDirectorySentry& operator=(const TH1AddDirectorySentry&);
  bool status_;
};


#endif
