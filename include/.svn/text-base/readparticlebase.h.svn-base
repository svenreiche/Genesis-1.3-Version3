#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>

#ifndef __GEN_READPARTICLEBASE__
#define __GEN_READPARTICLEBASE__

#include "genesis_fortran_common.h"
#include "genesis_fortran_calls.h"

using namespace std;

class ReadParticleBase{
 public:
  ReadParticleBase();
  virtual ~ReadParticleBase();
  virtual void init(string)=0;
  virtual void readpart(int)=0;
  //  virtual void finalize(int, bool)=0;
  //  virtual void writeRecord(int, bool)=0;
 protected:
  bool check();
  bool isOpen;
  double slicespacing,slicelength,refposition;
  int slicecount,beamletsize;
  double *work;
  int nwork;
};



#endif
