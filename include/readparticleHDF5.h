
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>


#ifndef __GEN_READPARTICLEHDF5__
#define __GEN_READPARTICLEHDF5__

#include "readparticlebase.h"
#include "hdf5.h"
#include "HDF5base.h"
#include "genesis_fortran_common.h"


using namespace std;

class ReadParticleHDF5 : public ReadParticleBase, public HDF5Base {
 public:
  ReadParticleHDF5();
  virtual ~ReadParticleHDF5();
  void init(string);
  void readpart(int islice);

  //  virtual void finalize(int, bool)=0;
  //  virtual void writeRecord(int, bool)=0;
 private:
  string filename;
  hid_t fid;
};



#endif

