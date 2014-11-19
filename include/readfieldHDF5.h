#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>


#ifndef __GEN_READFIELDHDF5__
#define __GEN_READFIELDHDF5__

#include "readfieldbase.h"
#include "hdf5.h"
#include "HDF5base.h"
#include "genesis_fortran_common.h"


using namespace std;

class ReadFieldHDF5 : public ReadFieldBase, public HDF5Base {
 public:
  ReadFieldHDF5();
  virtual ~ReadFieldHDF5();
  void init(string, int);
  void readfield(int islice);
  void readslippage();

  //  virtual void finalize(int, bool)=0;
  //  virtual void writeRecord(int, bool)=0;
 private:
  string filename;
  hid_t fid;
};



#endif

