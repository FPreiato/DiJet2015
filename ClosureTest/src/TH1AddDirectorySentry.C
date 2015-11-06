// system include files
#include "TH1.h"

// user include files
#include "TH1AddDirectorySentry.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TH1AddDirectorySentry::TH1AddDirectorySentry():
  status_(TH1::AddDirectoryStatus())
{
  TH1::AddDirectory(true);
}

TH1AddDirectorySentry::~TH1AddDirectorySentry()
{
  TH1::AddDirectory(status_);
}

