#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>


#ifndef __GEN_READPARTICLEBIN__
#define __GEN_READPARTICLEBIN__

#include "readparticlebase.h"
#include "genesis_fortran_common.h"


using namespace std;

class ReadParticleBinary : public ReadParticleBase{
 public:
  ReadParticleBinary();
  virtual ~ReadParticleBinary();
  void init(string);
  //  virtual void finalize(int, bool)=0;
  //  virtual void writeRecord(int, bool)=0;
 protected:
  string filename;
  ifstream fin;
  //  int firstout;
  //  bool dumppar,dumpfld;
};



#endif
